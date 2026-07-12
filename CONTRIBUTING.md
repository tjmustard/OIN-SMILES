# Contributing to OIN-SMILES

Hi there! We're thrilled that you'd like to contribute to this project. Contributions are released to the public under the [project's open source license](LICENSE).

Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

## Prerequisites for running and testing code

These are one-time installations required to test your changes locally as part of the pull request submission process.

1. Install [Python 3.10+](https://www.python.org/downloads/)
2. Install [uv](https://docs.astral.sh/uv/) for package management (`uv sync` installs all dependencies)

## Submitting a pull request

> [!NOTE]
> Please adhere to the project's primary objectives. If your pull request modifies core conversion or generation logic, it must pass all integration tests, including round-trip verification.

1. Fork and clone the repository.
2. Install the dependencies: `uv sync`.
3. Make your change.
4. Run the test suite: `uv run python -m unittest discover tests`.
5. Push to your fork and submit a pull request.

Here are a few things you can do that will increase the likelihood of your pull request being accepted:

- Follow the project's coding conventions (PEP 8; type hints on all public function signatures; docstrings for all public modules and functions; keep parsing logic separated from file I/O).
- Write tests for new functionality — unit tests in `tests/unit/`, integration tests in `tests/integration/`.
- Write a clear, descriptive commit message that explains the *why* rather than the *what*.
- For changes to the XYZ→OIN or OIN→XYZ pipelines, include a round-trip verification result showing RMSD < 1.0 Å on at least one curated example.

## AI contributions

> [!IMPORTANT]
> If you are using **any kind of AI assistance** to contribute to this project, it must be disclosed in the pull request or issue.

We welcome and encourage the use of AI tools (such as Claude, GitHub Copilot, Cursor, etc.) to help improve the codebase!

That being said, if you are using AI assistance, **this must be disclosed in the pull request**, along with the extent to which AI assistance was used.

An example disclosure:
> This PR was authored with Claude Code assistance. All changes were manually reviewed and verified.

Failure to disclose this makes it difficult to determine how much scrutiny to apply to the contribution.

When submitting AI-assisted contributions, please ensure they include:
- **Clear disclosure of AI use**.
- **Human verification** — you have manually verified that your changes work as intended.

## Resources

- [How to Contribute to Open Source](https://opensource.guide/how-to-contribute/)
- [Open Isomer Notation specification](README.md#oin-v36-inline-format)
