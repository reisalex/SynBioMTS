"""
Model interface to import models, write model wrappers, collect model output
and run statistics specific to model type

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

import stats


class interface():
    __name__ = "models.interface"

    def __init__(self):
        pass

    def RBS_Calculator_v1_0(self):
        pass

    def RBS_Calculator_v1_1(self):
        pass

    def RBS_Calculator_v2_0(self):
        pass

    def RBS_Calculator_v2_1(self):
        pass

    def UTR_Designer(self):
        pass

    def RBS_Designer(self):
        pass

    def EMOPEC(self):
        pass

    def calc_stats(self):
        pass