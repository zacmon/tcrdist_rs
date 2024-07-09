Developer tools
===============

Install `Rust <https://www.rust-lang.org/>`_ if necessary::

    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

Download the `maturin <https://github.com/PyO3/maturin>`_ library, assuming you are in an environment which has Python::

    pip install maturin

To compile tcrdist_rs, run::

    maturin develop --release --profile release
