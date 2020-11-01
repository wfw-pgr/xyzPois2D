import numpy as np


# ========================================================= #
# ===  make__sample.py                                  === #
# ========================================================= #

def make__sample():

    x_,y_,z_      = 0, 1, 2
    outFile       = "dat/source.dat"
    import nkUtilities.equiSpaceGrid as esg
    x1MinMaxNum   = [ 0.0, 1.0, 11 ]
    x2MinMaxNum   = [ 0.0, 1.0, 11 ]
    x3MinMaxNum   = [ 0.0, 0.0,  1 ]
    ret           = esg.equiSpaceGrid( x1MinMaxNum=x1MinMaxNum, x2MinMaxNum=x2MinMaxNum, \
                                       x3MinMaxNum=x3MinMaxNum, returnType = "structured" )
    aval, bval          = 1.0, 2.0
    ret[:,:,:,z_]       = aval * ret[:,:,:,x_] + bval * ret[:,:,:,y_]
    ret[:,1:-1,1:-1,z_] = 3.0 * ret[:,1:-1,1:-1,z_]**2
    ret[:, 0, :,z_]     = 0.0
    ret[:,-1, :,z_]     = 0.0
    ret[:, :, 0,z_]     = 0.0
    ret[:, :,-1,z_]     = 0.0
    
    import nkUtilities.save__pointFile as spf
    spf.save__pointFile( outFile=outFile, Data=ret )
    
    return()


# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #

if ( __name__=="__main__" ):
    make__sample()
