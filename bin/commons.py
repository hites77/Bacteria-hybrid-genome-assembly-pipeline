#!/usr/bin/env python3


def make_flag(name):
    return "--" + name


def bool_to_str(boolean):
    if not isinstance(boolean, bool):
        raise ValueError("boolean must be a bool")
    return "yes" if boolean else "no"
