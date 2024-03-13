1. Fork this repository
2. Create a branch for your fix of an issue/feature/doc (`git checkout -b my-fix`)
3. Stage your changes, run `pre-commit`, and run `pytest` until every thing looks fine
4. Commit your changes (`git commit -am 'Added some feature/fix/doc'`)
5. Push to the branch (`git push origin my-fix`)
6. Create new Pull Request

For code quality control, we use `pre-commit` hooks for formatting and type hints.

Try your best to add test cases for bug fixes and new features. We use `pytest` for testing.

For documentation, we should try to follow the [Google Python docstrings style](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html).
