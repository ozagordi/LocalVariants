#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    use_scm_version=True,
    setup_requires=['setuptools_scm', 'setuptools_scm_git_archive',
                    'isort', 'pytest-runner'],
    install_requires=['setuptools_scm', 'biopython'],
    tests_require=['pytest', 'flake8'],
    name='MinVar',
    description='A tool to analyse support files produced by `diri_sampler`',
    url='https://github.com/ozagordi/LocalVariants',
    author='Osvaldo Zagordi',
    author_email='firstname.lastname@gmail.com',
    packages=find_packages('src'),  # include all packages under src
    package_dir={'': 'src'},  # tell setuptools packages are under src
    package_data={'minvar': ['src/localvariants/db/*']},
    entry_points={
        'console_scripts': ['localvariants = localvariants.cli:main']
        },
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',

        # Pick your license as you wish (should match "license" above)
        'License :: Other/Proprietary License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        # 'Programming Language :: Python :: 2',
        # 'Programming Language :: Python :: 2.6',
        # 'Programming Language :: Python :: 2.7',
        # 'Programming Language :: Python :: 3',
        # 'Programming Language :: Python :: 3.2',
        # 'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.5']
)
