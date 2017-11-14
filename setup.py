from setuptools import setup, find_packages

setup(name='SnapGeneFileReader',
      version='0.1.9',
      author='yishaluo',
      author_email='yishaluo@gmail.com',
      description='A Python project to read and write Snapgene *.dna into dict, json, and biopython object.',
      long_description=open('README.md').read(),
      license='MIT',
      keywords="microplate biology parser converter report",
      packages= find_packages(),
      install_requires= ['biopython', 'xmltodict'])