import cell;
class configuration(object):

    class runParams;
        def __init__(self):
            self.SimLength = 100
            # TODO

    class simParams:
        def __init__(self):
            self.Arec = 1
            self.Brec = 2

    class CellParams:
        def __init__(self):
            self.Arec = 3
            self.Brec = 4

    def __init__(self):
        self.simconfig = self.runParams();
        self.enviorconfig = self.SimParams();


