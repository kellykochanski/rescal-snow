# How to contribute <a name="contribute">

Thank you for taking the time to contribute! We're excited to see your work.

[How can I contribute?](#contribute)
 - [Reporting bugs](#bugs)
 - [Adding features](#features)
    - [Git workflow and pull requests](#workflow)
 [Code of conduct](#conduct)
 
 ## Reporting bugs <a name="bugs">
 Please report bugs as Github Issues. Include steps necessary to reproduce the bug, any error text, and a description of the system on which the bug occurred. 
 For installation issues, include version numbers for the compiler, cmake, and the library dependencies.
 
 If you think you can fix the bug, great! See the features and pull request guidelines below for how to incorporate your bugfix into the Rescal-snow github.
 
 ## Adding features <a name="features">
 Rescal-snow grows by adding new features. Features include:
  - New capabilities for the Rescal-snow simulation (e.g. cohesion, or new inputs)
  - New workflow management, analysis, or visualization scripts

 The ideal feature:
  - Satisfies a clearly-defined purpose (e.g. "Add time-varying snowfall capability")
      - So "Fixed some typos in X.md" is a separate feature from "Added new tool for..."
  - Is self-contained in organized files, classes, and/or functions
  - Does not add dependencies to Rescal-snow or any existing script (optional dependencies OK)
  - Does not remove existing capabilities or break existing tutorials
      
We want your code to be be readable (so we can accept it), and described in the docs (so it will be found and used).
To make this possible, features that add new capabilities to Rescal-snow should include the following documentation:
  - Docstrings on every multi-line function
  - Author, date, and high-level purpose information in each new script/class
  - A tutorial or example in `docs` explaining the use of any novel capabilities
  - A reference that makes new capabilities discoverable through either `README.md` or `docs/rescal-snow-tutorial.md`
  - For particularly significant changes (project >1wk), appropriate additions to `NEWS.md` and `AUTHORS.md`
 
 ### Git workflow and pull requests <a name="workflow">
 We use the git branching workflow described [here](https://nvie.com/posts/a-successful-git-branching-model/) (if new to git, see link for useful commands).
 
 This consists of 4 main types of branches:
  - The `master` branch, which is changed only by new releases
  - A `develop` branch, which collects stable inter-release changes
  - An occasional `rescal-snow-release-X` branch which collects changes near new releases
  - *Feature* and *bugfix* branches for active development
  
When you add code to Rescal-snow, we encourage you to fork the code, 
then create a bug or feature branch titled `name/purpose`, e.g. `kk/snowfall`.
```bash
$ git checkout -b name/purpose develop
```
and commit your changes.

When the feature or bugfix is ready to be incorporated into Rescal-snow, open a [pull request](https://help.github.com/en/articles/creating-a-pull-request) against `develop`.
The feature should pass the tests in [test/test.sh](test/test.sh).
 
