import setuptools
import os

script_dir = os.path.join("src", "scripts")
# Collect entry scripts, annotated to keep the element type known.
# Use a quoted annotation so it is not evaluated at runtime (Python 3.6
# does not support subscripting the builtin `list`).
scripts: "list[str]" = [
    os.path.join(script_dir, fn)
    for fn in os.listdir(script_dir)
    if fn.endswith(".py")
]

with open("README.md", "r") as fh:
    long_description = fh.read()

# The returned Distribution object is intentionally discarded.
setuptools.setup(  # pyright: ignore[reportUnusedCallResult]
    name='pyBioInfo',
    version='1.0.1',
    author='Zonggui Chen',
    author_email='zgchen9@163.com',
    description='',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://gitee.com/Ckenen/pyBioInfo',
    package_dir={'': 'src'},
    packages=setuptools.find_packages(where="src"),
    package_data={'pyBioInfo': ['py.typed']},
    include_package_data=True,
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
