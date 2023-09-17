import setuptools
import os

script_dir = os.path.join("src", "scripts")
scripts = []
for fn in os.listdir(script_dir):
    if fn.endswith(".py"):
        scripts.append(os.path.join(script_dir, fn))

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='pyBioInfo',
    version='1.0.1',
    author='Zonggui Chen',
    author_email='ggchenzonggui@qq.com',
    description='',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://gitee.com/chenzonggui/pyBioInfo',
    package_dir={'': 'src'},
    packages=setuptools.find_packages(where="src"),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix"
    ],
    scripts=scripts,
    test_suite="tests",
    python_requires=">=3.6",
    install_requires=['biopython', 'pysam', 'numpy'],
)
