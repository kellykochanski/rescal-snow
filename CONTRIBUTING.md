#How to contribute <a name="contribute">

Thank you for taking the time to contribute! We're excited to see your work.

[How can I contribute?](#contribute)
 - [Reporting bugs](#bugs)
 - [Adding features](#features)
    - [Git workflow and pull requests](#workflow)
 [Code of conduct](#conduct)
 
 <a name="bugs">
 ##Reporting bugs
 Please report bugs as Github Issues. Include steps necessary to reproduce the bug, any error text, and a description of the system on which the bug occurred. 
 For installation issues, include version numbers for the compiler, cmake, and the library dependencies.
 
 If you think you can fix the bug, great! See the features and pull request guidelines below for how to incorporate your bugfix into the Rescal-snow github.
 
 <a name="features">
 ##Adding features
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
  
 <a name="workflow">
 ###Git workflow and pull requests
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
 
 <a name="conduct">
 #Contributor Covenant Code of Conduct 
##Our Pledge
In the interest of fostering an open and welcoming environment, we as contributors and maintainers pledge to make
participation in our project and our community a harassment-free experience for everyone, regardless of age, body size,
disability, ethnicity, sex characteristics, gender identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity and orientation. 

##Our Standards 
Examples of behavior that contributes to creating a positive environment include: 
- Using welcoming and inclusive language
- Being respectful of differing viewpoints and experiences
- Gracefully accepting constructive criticism
- Focusing on what is best for the community
- Showing empathy towards other community members 

Examples of unacceptable behavior by participants include: 
- The use of sexualized language or imagery and unwelcome sexual attention or advances
- Trolling, insulting/derogatory comments, and personal or political attacks
- Public or private harassment
- Publishing others' private information, such as a physical or electronic address, without explicit permission
- Other conduct which could reasonably be considered inappropriate in a professional setting 

##Our Responsibilities 
- Project maintainers are responsible for clarifying the standards of acceptable behavior and are expected to take appropriate
and fair corrective action in response to any instances of unacceptable behavior. 
- Project maintainers have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues,
and other contributions that are not aligned to this Code of Conduct, or to ban temporarily or permanently any contributor for
other behaviors that they deem inappropriate, threatening, offensive, or harmful. 

##Scope 
This Code of Conduct applies within all project spaces, and it also applies when an individual is representing the project or its
community in public spaces. Examples of representing a project or community include using an official project e-mail address,
posting via an official social media account, or acting as an appointed representative at an online or offline event.
Representation of a project may be further defined and clarified by project maintainers. 

##Enforcement 
Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by contacting the project team at
kelly.kochanski@gmail.com. All complaints will be reviewed and investigated and will result in a response that is deemed
necessary and appropriate to the circumstances. The project team is obligated to maintain confidentiality with regard to the
reporter of an incident. Further details of specific enforcement policies may be posted separately. 

Project maintainers who do not follow or enforce the Code of Conduct in good faith may face temporary or permanent
repercussions as determined by other members of the project’s leadership. 

##Attribution 
This Code of Conduct is adapted from the Contributor Covenant, version 1.4, available at https://www.contributor
covenant.org/version/1/4/code-of-conduct.html 
For answers to common questions about this code of conduct, see https://www.contributor-covenant.org/fa
 
 
 
