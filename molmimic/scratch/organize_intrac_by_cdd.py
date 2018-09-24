import binascii
import zlib

from pyasn1.codec.ber import decoder

import pandas as ps
import pymssql
from joblib import Parallel, delayed

def decode_residues(res):
    residues = []
    if res.startswith("0x"):
        res = res[2:]
    try:
        res = binascii.unhexlify(res)
    except:
        pass

    try:
        code, rest = decoder.decode(zlib.decompress(res, 16 + zlib.MAX_WBITS))
    except Exception as e:
        print e
        pass

    for i in xrange(len(code)):
        c = code[i]
        range_from, range_to, gi = tuple([c[j] for j in range(len(c))])
        for x in xrange(range_from, range_to + 1):
            residues.append(str(x))
    
    return ",".join(residues)

def get_interfaces_for_cdd(cdd_label, cdd_id):
    conn = pymssql.connect("IBIS_DEV", "anyone", "allowed", "Intrac")
    cursor = conn.cursor()
    sql = ""
    for table in xrange(1,87):
        db = int((table-1)/5)+1
        sql += "SELECT * FROM Intrac{db}.dbo.Inferred{tbl} i, Intrac{db}.dbo.InferredFace{tbl} f WHERE i.inf_id = f.inf_id AND i.nbr_superfam_id={cdd} AND i.nbr_type=3".format(tbl=table, db=db, cdd=int(cdd_id))
        if table < 86:
            sql += " UNION ALL "
    cursor.execute(sql) #This is fast, iteration over cursor is slow
    with open("inferred/{}.csv".format(cdd_label.replace("/", "").replace("'", "\'")), "a+") as cddf:
        for row in cursor:
            res = decode_residues(row[-1])
            print >> cddf, "\t".join(map(str, row[:-1]))+","+res
            
if __name__ == "__main__":
    CDD = pd.read_csv("StructDomSfam.csv", usecols=["sfam_id", "label"]).drop_duplicates().dropna()
    Parallel(n_jobs=64)(delayed(get_interfaces_for_cdd)(cdd["label"], cdd["sfam_id"]) for _, cdd in CDD.iterrows()) 