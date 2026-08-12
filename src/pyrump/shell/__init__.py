"""RUMP's interactive command shell.

A thin layer over the pyRUMP library: every command handler parses its
arguments, calls into ``pyrump.*``, and mutates :class:`~pyrump.shell.session.Session`
state. No physics lives here.
"""
