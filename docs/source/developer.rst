Developer tools
===============

Install Rust if necessary::

    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

Download the `maturin` library, assuming you are in an environment which has Python::

    pip install maturin

To compile tcrdist_rs, run::

    maturin develop --release --profile release
