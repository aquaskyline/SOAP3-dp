// Install zlib
// apt-get install zlib1g zlib1g-dev 
//
// Compile 
// gcc -lz testDrive.c libbam.a
//

#include "win32/zlib.h"
#include "sam.h"
#include "stdio.h"
#include "stdlib.h"

inline int SAMIUint8ConcatUint8(uint8_t * data, int * curSize,
                                uint8_t key) {
    data[(*curSize)++] = key;
}
inline int SAMIUint8ConcatUint32(uint8_t * data, int * curSize,
                                uint32_t key) {
    int i;
    int len = sizeof(uint32_t) / sizeof(uint8_t);
    uint8_t * key_8 = (uint8_t *) &key;
    for (i=0;i<len;i++) {
        data[(*curSize)++] = key_8[i];
    }
}
inline int SAMIUint8ConcatString(uint8_t * data, int * curSize,
                                char * key, int len) {
    int i;
    for (i=0;i<len;i++) {
        data[(*curSize)++] = key[i];
    }
}

int main () {
    bam_header_t sheader;
    
    sheader.n_targets = 2;
    char testSeq1[10] = "HELLOSEQ";
    char testSeq2[10] = "HELLOSE2";
    sheader.target_name = (char**) malloc(sizeof(char*)*2);
    sheader.target_name[0] = testSeq1;
    sheader.target_name[1] = testSeq2;
    sheader.target_len = (uint32_t*) malloc(sizeof(uint32_t)*2);
    sheader.target_len[0] = 10;
    sheader.target_len[1] = 10;
    
    sheader.hash = NULL;
    sheader.hash = NULL;
    
    sheader.l_text=0;
    sheader.text = (char*) malloc(sizeof(char)*5);
    sheader.text = "@TEST";
    
    samfile_t * sptr = samopen("ASD","wh",&sheader);
    
    bam1_t salignment;
    salignment.core.tid = 1;
    salignment.core.pos = 2456;
    salignment.core.bin = bam_reg2bin(0,0);
    salignment.core.qual = 255; //quality
    salignment.core.l_qname = 6; //lenght of the query name
    salignment.core.flag = 0;
    salignment.core.n_cigar = 1;
    salignment.core.l_qseq = 8; //length of the read
    
    salignment.core.mtid = 0;
    salignment.core.mpos = 1357;
    salignment.core.isize = 1;
        
    salignment.l_aux = 0;
    salignment.data_len = 0;
    salignment.m_data = 255;
    salignment.data = malloc(sizeof(uint8_t)*255);
    
    //Name
    SAMIUint8ConcatString(salignment.data,&(salignment.data_len),"QNAME",6);
    //CIGAR
    SAMIUint8ConcatUint32(salignment.data,&(salignment.data_len),128);
    //Read
    SAMIUint8ConcatUint8(salignment.data,&(salignment.data_len),17);
    SAMIUint8ConcatUint8(salignment.data,&(salignment.data_len),17);
    SAMIUint8ConcatUint8(salignment.data,&(salignment.data_len),34);
    SAMIUint8ConcatUint8(salignment.data,&(salignment.data_len),34);
    //Quality
    SAMIUint8ConcatUint8(salignment.data,&(salignment.data_len),0xff);
    SAMIUint8ConcatUint8(salignment.data,&(salignment.data_len),0xff);
    SAMIUint8ConcatUint8(salignment.data,&(salignment.data_len),0xff);
    SAMIUint8ConcatUint8(salignment.data,&(salignment.data_len),0xff);
    SAMIUint8ConcatUint8(salignment.data,&(salignment.data_len),0xff);
    SAMIUint8ConcatUint8(salignment.data,&(salignment.data_len),0xff);
    SAMIUint8ConcatUint8(salignment.data,&(salignment.data_len),0xff);
    SAMIUint8ConcatUint8(salignment.data,&(salignment.data_len),0xff);
    
    unsigned int auxStart = salignment.data_len;
    SAMIUint8ConcatString(salignment.data,&(salignment.data_len),"XAZTESTING",11);
    salignment.l_aux += salignment.data_len - auxStart; 
    //salignment.data = "QNAME\0\0\0\0\0AAAACCCC!!!!!!!!ASD";
    
    
    unsigned int * asd = bam1_cigar(&salignment);
    
    printf("%u\n",(uint32_t) (*asd));
    
    samwrite(sptr,&salignment);
    
    salignment.core.tid = 0;
    salignment.core.pos = 1357;
    samwrite(sptr,&salignment);
    
    
    samclose(sptr);

    return 0;
}
