from setuptools import setup

setup(
    name='pySPHERE',
    version='0.1',

    description='VLT/SPHERE reduction and analysis code',
    url='https://github.com/avigan/pySPHERE',
    author='Arthur Vigan',
    author_email='arthur.vigan@lam.fr',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Professional Astronomers',
        'Topic :: High-contrast Imaging and Spectroscopy',
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: BSD License'
    ],
    keywords='vlt sphere ifs irdis reduction',
    packages=['pysphere', 'pysphere.utils'],
    install_requires=[
        'numpy', 'scipy', 'astropy', 'pandas', 'matplotlib'
    ],
    include_package_data=True,
    package_data={
        'pysphere': ['data/*.txt', 'data/*.dat', 'data/*.fits'],
    },
    zip_safe=False
)
