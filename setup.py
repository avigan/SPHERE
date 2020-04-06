from setuptools import setup

setup(
    name='pysphere',
    version='1.1',

    description='Reduction and analysis code for the VLT/SPHERE instrument',
    url='https://github.com/avigan/SPHERE',
    author='Arthur Vigan',
    author_email='arthur.vigan@lam.fr',
    license='MIT',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Professional Astronomers',
        'Topic :: High-contrast Imaging and Spectroscopy',
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License'
    ],
    keywords='vlt sphere ifs irdis reduction exoplanet pipeline',
    packages=['pyshere', 'pysphere.utils'],
    install_requires=[
        'numpy', 'scipy', 'astropy', 'pandas', 'matplotlib'
    ],
    include_package_data=True,
    package_data={
        'pysphere': ['data/*.txt', 'data/*.dat', 'data/*.fits',
                     'instruments/*.ini', 'instruments/*.dat'],
    },
    zip_safe=False
)
