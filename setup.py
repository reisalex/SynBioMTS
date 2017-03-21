"""
Python-packaging for testsfm

Copyright 2017 Alexander C. Reis, Howard M. Salis, all rights reserved.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

Please cite:

  Alexander C. Reis, and Howard M. Salis
  An automated model test system for systematic development and improvement of
  gene expression models, Nature Methods (2017)

"""

from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='testsfm',
      version='1.0',
      description='Test suite for DNA sequence-function models',
      url='http://github.com/reisalex/test-sfm',
      author='Alexander C. Reis',
      author_email='alex.reis@psu.edu',
      license='MIT',
      packages=['testsfm'],
      install_requires=[
          'numpy',
          'scipy',
          'pandas',
          'biopython'
      ],
      zip_safe=False)