import sys
import numpy                      as np
import nkUtilities.load__config   as lcf
import nkUtilities.cMapTri        as cmt
import nkUtilities.configSettings as cfs


# ========================================================= #
# ===  display                                          === #
# ========================================================= #
def display():
    # ------------------------------------------------- #
    # --- [1] Arguments                             --- #
    # ------------------------------------------------- #
    config   = lcf.load__config()
    datFile1 = "dat/result.dat"
    datFile2 = "dat/source.dat"
    pngFile1 = datFile1.replace( "dat", "png" )
    pngFile2 = datFile2.replace( "dat", "png" )

    # ------------------------------------------------- #
    # --- [2] Fetch Data                            --- #
    # ------------------------------------------------- #
    import nkUtilities.load__pointFile as lpf
    Data1  = lpf.load__pointFile( inpFile=datFile1, returnType="point" )
    xAxis1 = Data1[:,0]
    yAxis1 = Data1[:,1]
    zAxis1 = Data1[:,2]
    
    Data2  = lpf.load__pointFile( inpFile=datFile2, returnType="point" )
    xAxis2 = Data2[:,0]
    yAxis2 = Data2[:,1]
    zAxis2 = Data2[:,2]
    
    # ------------------------------------------------- #
    # --- [3] config Settings                       --- #
    # ------------------------------------------------- #
    cfs.configSettings( configType="cMap_def", config=config )
    config["FigSize"]        = (5,5)
    config["cmp_position"]   = [0.16,0.12,0.97,0.88]
    config["xTitle"]         = "X (m)"
    config["yTitle"]         = "Y (m)"
    config["cmp_xAutoRange"] = True
    config["cmp_yAutoRange"] = True
    config["cmp_xRange"]     = [-5.0,+5.0]
    config["cmp_yRange"]     = [-5.0,+5.0]

    # ------------------------------------------------- #
    # --- [4] plot Figure                           --- #
    # ------------------------------------------------- #
    cmt.cMapTri( xAxis=xAxis1, yAxis=yAxis1, cMap=zAxis1, pngFile=pngFile1, config=config )
    cmt.cMapTri( xAxis=xAxis2, yAxis=yAxis2, cMap=zAxis2, pngFile=pngFile2, config=config )


# ======================================== #
# ===  実行部                          === #
# ======================================== #
if ( __name__=="__main__" ):
    display()

