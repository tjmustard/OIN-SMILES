# Contributing to {{PROJECT_NAME}}

Hi there! We're thrilled that you'd like to contribute to {{PROJECT_NAME}}. Contributions to this project are released to the public under the [project's open source license](LICENSE).

Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

## Prerequisites for running and testing code

These are one-time installations required to be able to test your changes locally as part of the pull request (PR) submission process.

1. Install [Python 3.10+](https://www.python.org/downloads/)
2. Install [uv](https://docs.astral.sh/uv/) for package management

## Submitting a pull request

> [!NOTE]
> Please adhere to the project's primary objectives. If your pull request modifies core components, it must pass all integration tests.

1. Fork and clone the repository.
2. Install the dependencies: `uv sync`.
3. Make your change.
4. Run the integration test suite: `uv run pytest tests/`.
5. Push to your fork and submit a pull request.

Here are a few things you can do that will increase the likelihood of your pull request being accepted:

- Follow the project's coding conventions (type hints required for all functions, Google Style docstrings).
- Write tests for new functionality.
- Update the documentation suite using the `/hyper-document` workflow if you are utilizing the agent framework.
- Write a clear, descriptive commit message.

## AI contributions

> [!IMPORTANT]
> If you are using **any kind of AI assistance** to contribute to this project, it must be disclosed in the pull request or issue.

We welcome and encourage the use of AI tools (such as Claude, GitHub Copilot, Cursor, etc.) to help improve the codebase! 

That being said, if you are using AI assistance, **this must be disclosed in the pull request**, along with the extent to which AI assistance was used.

An example disclosure:
> This PR was authored by Claude using the Hypergraph Coding Agent Framework.

Failure to disclose this makes it difficult to determine how much scrutiny to apply to the contribution.

### What we're looking for

When submitting AI-assisted contributions, please ensure they include:
- **Clear disclosure of AI use**.
- **Human verification** - You have manually verified that your changes work as intended.

## Resources

- [Tutorial](./docs/TUTORIAL.md)
- [How to Contribute to Open Source](https://opensource.guide/how-to-contribute/)
