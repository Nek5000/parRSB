#include "gencon-impl.h"
#include "exa-memory.h"

int setGlobalID(exaHandle h,Mesh mesh){
  exaComm c  =exaGetComm(h);
  exaInt rank=exaCommRank(c);
  exaInt size=exaCommSize(c);

  exaInt nPoints=exaArrayGetSize(mesh->elements);
  Point   points=exaArrayGetPointer(mesh->elements);

  exaInt bin=1;
  if(nPoints==0) bin=0;
  exaComm nonZeroRanks;
  exaCommDup(&nonZeroRanks,c);
  exaCommSplit(&nonZeroRanks,bin);

  if(bin==1){
    exaLong count=0;
    exaInt i;
    for(i=0;i<nPoints;i++)
      if(points[i].ifSegment) count++;

    //for(int i=0;i<size;i++){
    //  exaCommBarrier(nonZeroRanks);
    //  if(i==rank)
    //    for(int j=0;j<nPoints;j++)
    //      printf("rank=%d segment[%02d]=%d\n",i,j,points[j].ifSegment);
    //}

    exaLong out[2][1],buff[2][1],in[1];
    in[0]=count+!rank;
    exaCommScan(nonZeroRanks,out,in,buff,1,exaLong_t,exaAddOp);
    exaLong start=out[0][0];
    exaDebug(h,"rank=%d start=%lld size=%d\n",rank,start,in[0]);

    start-=(rank>0?1:0);
    count=0;
    for(i=0;i<nPoints;i++){
      if(points[i].ifSegment) count++;
      points[i].globalId=start+count;
    }
  }

  exaDestroy(nonZeroRanks);
  return 0;
}

int sendBack(exaHandle h,Mesh mesh){
  exaComm c  =exaGetComm(h);
  exaInt rank=exaCommRank(c);
  exaInt size=exaCommSize(c);

  exaInt nPoints=exaArrayGetSize(mesh->elements);
  exaArrayTransfer(mesh->elements,offsetof(struct Point_private,origin),1,c);
  exaSortArray(mesh->elements,exaLong_t,offsetof(struct Point_private,sequenceId));
}
