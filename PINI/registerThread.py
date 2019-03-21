# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 14:10:35 2016

@author: broche
"""

# -*- coding: utf-8 -*-

import PyMcaQt as qt

class registerThread(qt.QThread):
    RegDone = qt.pyqtSignal()

    def __init__(self,register,parent):

        qt.QThread.__init__(self, parent)
        self.daddy=parent
        self.register=register

    def run(self):
        self.register.Execute()
        self.RegDone.emit()