class InsaneBlockOnTheLooseException( Exception ):
    def __init__(self, badBlock):
        Exception.__init__(self, "Error, an insane block on the loose! Call animal control, and tell them to bring a shotgun and a net, and have them search the following, offending block!\n Insane Block:\n\t " + str(badBlock))
        self.badBlock = badBlock[:]
        
    
