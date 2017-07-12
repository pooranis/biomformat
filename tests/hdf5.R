library(rhdf5)
setwd("~/git/biomformat/")
file <- "testdata/rich_dense_otu_table.v2.biom"

loc <- H5Fopen(file)
ga <- H5Aopen(loc, "type")
ar <- H5Aread(ga)
print(ar)
stop()

loc <- H5Fopen(file)
h5space4 = H5Screate("H5S_SCALAR")
#h5space4 = H5Screate_simple(0,NULL)
typeid <- H5Tcopy("H5T_C_S1")
H5Tset_size(typeid, size=10)
ga <- H5Acreate(loc, "type", typeid, h5space4)
ccc <- H5Awrite(ga, "OTU table")
H5Sclose(h5space4)
H5Aclose(ga)
H5close()
stop()

loc <- H5Fopen(file)
H5Adelete(loc, "type")
H5close()

# create a file and write something
h5createFile("ex_H5A.h5")
h5write(1:15, "ex_H5A.h5","A")

# write an attribute 'unit' to 'A'
fid <- H5Fopen("ex_H5A.h5")
did <- H5Dopen(fid, "A")
sid <- H5Screate_simple(c(1,1))
tid <- H5Tcopy("H5T_C_S1")

H5Tset_size(tid, 10L)
aid <- H5Acreate(did, "unit", tid, sid)
aid
H5Awrite(aid, "liter")
H5Aclose(aid)
H5Sclose(sid)
H5Aexists(did, "unit")
H5Dclose(did)
H5Fclose(fid)
h5dump("ex_H5A.h5")