#this function takes a single input, which is a table with three colums. Column 1 contains the chromosome (with no chr prefix),
#column 2 contain the last ukb WGS block for each chromosome (e.g. b4979 for chromosome 1) Column 3 contains the last base number 
#for each chromosome (i.e. the chromosome length). WGS blocks are numbered in a 0-based manner (i.e. chromosome 1 has 4980 blocks). 
#Each block contains 50kb of sequence, with the exception of regions 49100001-49150000 of chromosome 4 and 41850001-41900000 of 
#chromosome 10, which have 5kb blocks. The function creates a table which functions as map of exact coordinates contained by 
#each WGS block.

create_wgs_block_map <- function(block_index_filepath){
  
  index = read.table(block_index_filepath,header=T,sep="\t")
  index[,"LAST_BLOCK"] = gsub("b","",index[,"LAST_BLOCK"])
  index[,"LAST_BLOCK"] = as.numeric(index[,"LAST_BLOCK"])
  
  start = c()
  end = c()
  chr = c()
  wgs_block = c()
  file_name = c()
  
  for (c in c(1:22,"X")){
    print(paste("Processing chromosome ",c,"...",sep=""))
    
    nblocks = index[which(index["CHROM"]==c),"LAST_BLOCK"] #0-based
    clength = index[which(index["CHROM"]==c),"LAST_BASE"]
    
    chr = c(chr,rep(c,nblocks+1)) #the +1 is needed because blocks are numbered in a 0-based manner
    for(b in 0:nblocks){
      
      #define start and end coordinates for each block
      s = (b*50000)+1
      e = (b+1)*50000
      
      if ((c==4 & b==982) | (c==10 & b==837)){ #first small blocks on c4 and c10 got a normal 50kb-based start and a 5kb-based end, which gets redefined
        e = (b*50000)+5000
      }
      else if (c==4 & b>982 & b<992){ #blocks between 983 and 991 have a 5kb-based start and a 5kb-based end - both get redefined
        s = ((982*50000)+5000*(b-982))+1
        e = s+4999
      }
      else if (c==10 & b>837 & b<847){ #same for blocks between 838 and 846 on chr10
        s = ((837*50000)+5000*(b-837))+1
        e = s+4999
      }
      else if ((c==4 & b>=992) | (c==10 & b>=847)){ #blocks after the small ones are normally 50kb based, but with shifted coordinates with respect to the b number
        s = ((b-9)*50000)+1
        e = (b+1-9)*50000
      }
      else{
        #do nothing
      }
      
      if(e >= clength){
        e = clength
      }
      else{
        #do nothing
      }
      
      #define block name for each block
      wb = paste("b",b,sep="")
      
      #define file name for each block
      fn = paste("ukb24304_c",c,"_",wb,"_v1.vcf.gz",sep="")
      
      start = c(start,s)
      end = c(end,e)
      wgs_block = c(wgs_block,wb)
      file_name = c(file_name,fn)
    }
  }
  block_map = cbind(chr,start,end,wgs_block,file_name)
  write.table(block_map,"WGS_200k_block_map.tsv",col.names=T,row.names=F,sep="\t")
  }