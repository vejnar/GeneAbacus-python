# -*- coding: utf-8 -*-

#
# Copyright (C) 2015-2021 Charles E. Vejnar
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
#

import json
import os
import zlib

import lz4.frame
import numpy as np

def get_flength(features):
    return sum([e-s for s, e in features])

def parse_feat_line(line):
    cells = line.rstrip().split('\t')
    return [cells[0], int(cells[1])]

def pfopen(filename, path_features=None, fon_name='transcript_stable_id', fon_coords='exons', features=[], rw=False):
    assert os.path.exists(filename), f'{filename} not found'
    assert os.path.exists(path_features), f'{path_features} not found'
    # Open file
    if filename.endswith('.bin'):
        f = open(filename, 'rb')
    elif filename.endswith('.bin.lz4'):
        f = lz4.frame.open(filename, 'rb')
    else:
        raise ValueError('Invalid format')
    # Open feature file
    if path_features.endswith('.fon1.json'):
        features = [[ft[fon_name], get_flength(ft[fon_coords])] for ft in json.load(open(path_features, 'rt'))['features']]
    elif path_features.endswith('.tab'):
        with open(path_features, 'rt') as ftab:
            features = [parse_feat_line(l) for l in ftab]
    else:
        raise ValueError(f'Unknown format: {path_features}')
    # Read version
    detected_version = np.frombuffer(f.read(1), dtype='uint8')[0]
    assert detected_version in [3], f'Unknown version: {detected_version}'
    # Read Length
    full_length = np.frombuffer(f.read(4), dtype='uint32')[0]
    # Read checksum
    detected_csum = np.frombuffer(f.read(4), dtype='uint32')[0]
    # Compute Adler-32 checksum
    computed_csum = zlib.adler32(np.array([ft[1] for ft in features], dtype='uint32').tobytes())
    # Check checksum
    assert detected_csum == computed_csum, f'Feature file does not correspond to profile checksum ({detected_csum} != {computed_csum})'
    # Read file into Numpy array
    raw_profiles = np.zeros(full_length, dtype='float32')
    f.readinto(raw_profiles.data)
    if rw == False:
        raw_profiles.flags.writeable = False
    f.close()
    # Split array per feature
    profiles = {}
    i = 0
    for name, length in features:
        profiles[name] = raw_profiles[i:i+length]
        i += length
    assert i == len(raw_profiles), f'{i} out of {len(raw_profiles)} was read from input'
    return profiles

def pfwrite(profiles, filename, path_features=None, fon_name='transcript_stable_id', fon_coords='exons', features=[], pformat='binary'):
    # Open feature
    if path_features is not None:
        assert os.path.exists(path_features), f'{path_features} not found'
        if path_features.endswith('.fon1.json'):
            features = [[ft[fon_name], get_flength(ft[fon_coords])] for ft in json.load(open(path_features, 'rt'))['features']]
        elif path_features.endswith('.tab'):
            with open(path_features, 'rt') as ftab:
                features = [parse_feat_line(l) for l in ftab]
        else:
            raise ValueError(f'Unknown format: {path_features}')
    elif len(features) == 0:
        raise ValueError('No input features')
    if pformat == 'binary':
        # Open file
        if filename.endswith('.bin'):
            f = open(filename, 'wb')
        elif filename.endswith('.bin.lz4'):
            f = lz4.frame.open(filename, 'wb')
        else:
            raise ValueError('Invalid format')
        # Write current version
        f.write(np.uint8(3).tobytes())
        # Write length
        f.write(np.uint32(sum([ft[1] for ft in features])).tobytes())
        # Compute Adler-32 checksum
        computed_csum = zlib.adler32(np.array([ft[1] for ft in features], dtype='uint32').tobytes())
        # Write checksum
        f.write(np.uint32(computed_csum).tobytes())
        # Write profiles
        for ft in features:
            p = profiles[ft[0]]
            if len(p) == ft[1]:
                if p.dtype != 'float32':
                    q = p.astype('float32')
                else:
                    q = p
                f.write(q.tobytes())
            else:
                raise ValueError(f'Invalid length: {ft[0]} ({len(p)} != {ft[1]})')
        f.close()
    elif pformat == 'csv':
        # Open file
        if filename.endswith('.csv'):
            f = open(filename, 'wt')
        elif filename.endswith('.csv.lz4'):
            f = lz4.frame.open(filename, 'wt')
        else:
            raise ValueError('Invalid format')
        # Write profiles
        for ft in features:
            f.write(','.join([ft[0], ft[1], ' '.join([f'{i:.5f}' for i in profiles[ft[0]]])]))
            f.write('\n')
        f.close()
