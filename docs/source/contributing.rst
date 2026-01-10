============
Contributing
============

We welcome contributions to HPC-ScaleTest!

How to Contribute
=================

Code Contributions
------------------

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests
5. Submit a pull request

Documentation
-------------

Documentation improvements are always welcome:

* Fix typos
* Add examples
* Clarify explanations
* Add tutorials

Bug Reports
-----------

Report bugs via GitHub Issues:

* Describe the problem
* Include steps to reproduce
* Provide system information
* Attach relevant logs

Feature Requests
----------------

Suggest new features via GitHub Issues:

* Describe the use case
* Explain expected behavior
* Consider implementation approach

Development Workflow
====================

Setup Development Environment
------------------------------

.. code-block:: bash

   git clone https://github.com/user/hpc-scaletest.git
   cd hpc-scaletest
   python3 -m venv venv
   source venv/bin/activate
   pip install -e ".[dev]"

Create Feature Branch
---------------------

.. code-block:: bash

   git checkout -b feature/my-new-feature

Make Changes
------------

Follow code style guidelines (see :doc:`developer_guide`).

Run Tests
---------

.. code-block:: bash

   pytest tests/ -v
   black .
   flake8 .

Commit Changes
--------------

Write clear commit messages:

.. code-block:: bash

   git commit -m "Add feature X for Y"

Push and Create PR
------------------

.. code-block:: bash

   git push origin feature/my-new-feature

Then create a pull request on GitHub.

Coding Standards
================

* Follow PEP 8
* Use type hints
* Write docstrings
* Add unit tests
* Update documentation

Pull Request Guidelines
=======================

* One feature per PR
* Include tests
* Update documentation
* Pass CI checks
* Respond to review comments

Community Guidelines
====================

* Be respectful
* Be constructive
* Be patient
* Help others

License
=======

By contributing, you agree that your contributions will be licensed under the MIT License.

See Also
========

* :doc:`developer_guide` - Developer documentation
* :doc:`architecture` - Architecture overview
