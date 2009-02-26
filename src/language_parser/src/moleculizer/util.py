# These are things that shoud be defined in Main.

def And( x, y):
    return x and y



def DataUnifier( anArray):
    # This function takes an array that is a bunch of comma ended and semicolon ended 
    # lines, and combines lines so that there is one statement per line.

    # ['GDP,', 'mass = 100.0;', 'GTP,', 'mass = 110.0;', 'Ste2(to-Ste4,  to-alpha),', 'mass = 100.0,', 'preferred_state = top, front;', 'Ste4(to-Ste2, to-Gpa1, to-Ste5 {default, obstructed}, to-Fus3, *ModSite {None, Phosphorylated, DoublyPhosphorylated} ),', 'mass = 100.0;', 'Ste5( to-Ste4, to-Ste11, to-Ste7, to-Fus3),', 'mass = 100.0;', 'Gpa1( to-Ste4 { default, GTP-bound-shape}, to-GXP { default, bound-shape} ),', 'mass = 100.0;', 'alpha,', 'mass = 100.0;']
    if (len(anArray) <= 1): return anArray

    newArray = []
    current_string = ""

    for ndx in range(len(anArray)):
        current_string += anArray[ndx].strip()
        current_string = current_string.strip()
        if current_string[-1] == ";":
            newArray.append(current_string)
            current_string = ""

    return newArray

