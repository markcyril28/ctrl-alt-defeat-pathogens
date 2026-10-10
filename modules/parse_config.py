#!/usr/bin/env python3
"""Config file utilities for the Ologist Workshop pipeline.

Subcommands:
    to-json <input.toml> <output.json>   Convert TOML config to JSON
    species-name <config.json>           Print the species name from JSON config
    species-dir <config.json>            Print the species directory from JSON config
    accessions <config.json>             Print space-separated fetch target accessions
    docking-recipes <input.toml>         Print TSV docking recipes from TOML config
    species-name-toml <input.toml>       Print species name directly from TOML config
"""
import sys
import json


def _load_toml(path):
    try:
        import tomllib
    except ImportError:
        import tomli as tomllib
    with open(path, 'rb') as f:
        return tomllib.load(f)


def cmd_to_json(args):
    cfg = _load_toml(args[0])
    with open(args[1], 'w') as f:
        json.dump(cfg, f)


def cmd_species_name(args):
    c = json.load(open(args[0]))
    print(c['species']['name'])


def cmd_species_dir(args):
    c = json.load(open(args[0]))
    print(c['species']['dir'])


def cmd_accessions(args):
    c = json.load(open(args[0]))
    print(' '.join(t['accession'] for t in c.get('fetch_target', [])))


def cmd_docking_recipes(args):
    cfg = _load_toml(args[0])
    sd = cfg['species']['dir']
    for t in cfg.get('dock_target', []):
        ps = t.get('pre_selection', '-')
        sr = t.get('split_resn', '-')
        km = 'yes' if t.get('keep_metals', False) else 'no'
        ctr = t.get('center', 'user')
        box = t.get('box', 22)
        rr = t.get('redock_resn', '-')
        note = t.get('note', '')
        rec = '{}/{}'.format(sd, t['receptor'])
        print('\t'.join([t['key'], rec, ps, sr, km, str(ctr), str(box), rr, note]))


def cmd_species_name_toml(args):
    cfg = _load_toml(args[0])
    print(cfg['species']['name'])


if __name__ == '__main__':
    commands = {
        'to-json': cmd_to_json,
        'species-name': cmd_species_name,
        'species-dir': cmd_species_dir,
        'accessions': cmd_accessions,
        'docking-recipes': cmd_docking_recipes,
        'species-name-toml': cmd_species_name_toml,
    }
    if len(sys.argv) < 2 or sys.argv[1] not in commands:
        print("usage: parse_config.py <command> [args...]", file=sys.stderr)
        print("commands: " + ", ".join(commands), file=sys.stderr)
        sys.exit(1)
    commands[sys.argv[1]](sys.argv[2:])
