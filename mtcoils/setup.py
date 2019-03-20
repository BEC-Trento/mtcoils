#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Create: 03-2019 - Carmelo Mordini <carmelo> <carmelo.mordini@unitn.it>

"""
Predefined setup loading
"""
from pathlib import Path
from ruamel.yaml import YAML
yaml = YAML(typ='safe')

data = Path(__file__).parent / 'data' / 'magnetic-trap-bec1.yaml'

with open(data) as f:
    setup_bec1 = yaml.load(f)
