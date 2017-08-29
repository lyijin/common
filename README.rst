==================================================
"common": scripts commonly used across my projects
==================================================

Personal preferences and coding styles are listed below.

Python
------
- Python 3. Some scripts use features that was made possible in later versions (e.g. print (\*foo, \*bar), only possible in Python 3.4). For best results, use the *latest* version of Python 3.
- Follows PEP8 as closely as possible, 80 chars per line, with some exceptions.
- Four spaces, not tabs.
- Imports: libraries in standard library (e.g. ``argparse``, ``collections``, ``csv``, ...) first, followed by very common scientific Python libraries (e.g. ``numpy``, ``scipy``, ...), then my own custom scripts (e.g. ``natural_sort``, ``parse_*``, ...).

Shell scripts
-------------
- Targeted towards bash. No idea whether they'd work in zsh.

Input files
-----------
Unless otherwise specified...

- Plain-text files that are tab-delimited with UNIX line endings ("\\n").
- ASCII-encoded; though UTF-8 should work with no modifications because Python 3 is awesome.

Additional notes
----------------
While it might appear disorganised, the prefix of a filename usually gives one a good idea of what the script is capable of doing. For instance:
- ``parse_*.py``: parses ``*``, returns stuff. Usually written in a way that allows one to import that script as a library in a Python script.
- ``plot_*.py``: plots ``*``. Not publication-quality, it's more to visualise stuff that Excel cannot handle (e.g. files with millions of lines).

My scripts were written over the course of many many years, and I try to make it as bug-free as I could--but I'm pretty sure they're still some edge cases that I've never envisioned. Oh well.

Also, as I got more experienced I started using shortcuts e.g. list comprehensions and splatting, so my code might be hard to understand at times. I have however tried to comment my code to aid understanding--if there's something I learnt while coding in Perl, it's that if I didn't comment religiously in my code, I'd forget what I've done in six months. Or perhaps I'm just terrible at Perl!
