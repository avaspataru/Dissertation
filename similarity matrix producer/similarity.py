
import numpy as np
import pandas as pd
from sklearn.preprocessing import scale

def readData(fileName):
    fp = open( fileName, 'r')

    line = fp.readline() #ignore headers
    cnt = 1
    array = []
    while line:
       line = fp.readline()
       if line == "":
           continue
       qcnt = 0 #number of "
       dcnt = 0 #number of $
       cluster_id_s = ""
       gene = ""
       avg_expr_s = ""
       for c in line:
          if c == '"':
              qcnt=qcnt+1
              continue
          if c == '$':
              dcnt=dcnt+1
              continue
          if c== ' ':
              continue
          if(qcnt == 3 and dcnt == 0): #in gene name
              gene = gene +c
          if(qcnt == 3 and dcnt == 1 ): #in cluster id
              cluster_id_s = cluster_id_s + c
          if(qcnt == 3 and dcnt == 2 ): # in avg value
              avg_expr_s = avg_expr_s + c

       avg_expr = float(avg_expr_s)
       cluster_id = int(cluster_id_s)

       elem = [gene, cluster_id, avg_expr]
       array.append(elem)
       cnt += 1
    fp.close()

    d = ['gene', 'cluster_id', 'avg_expr']
    df = pd.DataFrame(array, columns=d)
    return df

def main():
    print "Computing similarities for pre and post identified clusters."
    pre_data = readData('preGenesClusters.txt')
    post_data = readData('postGenesClusters.txt')
    pre_clusters = pre_data['cluster_id'].unique()
    post_clusters = post_data['cluster_id'].unique()

    top_n = 10
    if(top_n ==-1):
        for i in pre_clusters:
            print "-------------------------------------------------------------------"
            losses = []
            for j in post_clusters:
                pre_cluster = pre_data.loc[pre_data['cluster_id'] == i]
                post_cluster = post_data.loc[post_data['cluster_id']==j]

                lookup_genes = pre_cluster['gene'].unique()
                loss = 0
                ngene = lookup_genes.size
                for gene in lookup_genes:
                    avg_expr_pre = pre_cluster.loc[pre_cluster['gene'] == gene]['avg_expr'].item()
                    if(post_cluster.loc[post_cluster['gene'] == gene].empty): #gene doesn't exist in post cluster
                        loss = loss + avg_expr_pre*avg_expr_pre
                        continue
                    avg_expr_post = post_cluster.loc[post_cluster['gene'] == gene]['avg_expr'].item()
                    #print avg_expr_pre
                    loss = loss + (avg_expr_pre - avg_expr_post) *  (avg_expr_pre - avg_expr_post)

                lookup_genes = post_cluster['gene'].unique()
                for gene in lookup_genes:
                    avg_expr_pre = post_cluster.loc[post_cluster['gene'] == gene]['avg_expr'].item()
                    if(pre_cluster.loc[pre_cluster['gene'] == gene].empty): #gene doesn't exist in post cluster
                        loss = loss + avg_expr_post*avg_expr_post
                        ngene = ngene+ 1
                loss = loss / ngene
                losses = losses + [loss]
                print "compare pre_" + str(i) + ", post_" + str(j) + ": loss " + str(loss) + "; looked at " + str(ngene) + " genes."
            #losses = scale( losses, axis=0, with_mean=True, with_std=True, copy=True )
            #print losses
    else:
        for i in pre_clusters:
            print "-------------------------------------------------------------------"
            for j in post_clusters:
                pre_cluster = pre_data.loc[pre_data['cluster_id'] == i]
                post_cluster = post_data.loc[post_data['cluster_id']==j]

                top_genes = pre_cluster.nlargest(top_n,'avg_expr')
                lookup_genes = top_genes['gene'].unique()
                loss = 0
                ngene = lookup_genes.size
                for gene in lookup_genes:
                    avg_expr_pre = pre_cluster.loc[pre_cluster['gene'] == gene]['avg_expr'].item()
                    if(post_cluster.loc[post_cluster['gene'] == gene].empty): #gene doesn't exist in post cluster
                        loss = loss + avg_expr_pre*avg_expr_pre
                        continue
                    avg_expr_post = post_cluster.loc[post_cluster['gene'] == gene]['avg_expr'].item()

                    loss = loss + (avg_expr_pre - avg_expr_post) *  (avg_expr_pre - avg_expr_post)
                loss = loss / ngene
                print "compare pre_" + str(i) + ", post_" + str(j) + ": loss " + str(loss) + "; looked at " + str(ngene) + " genes."


main()
