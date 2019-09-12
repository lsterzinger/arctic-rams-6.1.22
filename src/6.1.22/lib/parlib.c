/*
!
! Copyright (C) 1991-2004  ; All Rights Reserved ; Colorado State University
! Colorado State University Research Foundation ; ATMET, LLC
!
! This file is free software; you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later version.
!
! This software is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with this
! code; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!======================================================================================
*/

#include <stdio.h>
#include "utils_sub_names.h"

#if defined (RAMS_MPI)
  #include <mpi.h>
  #define MAX_MSGTAG 1000
#else
  #define MAX_MSGTAG 1
  #define MPI_Request int
#endif

int flag_msgtag=0;
MPI_Request mpi_msgtags[MAX_MSGTAG];

char *ibuff;
int ipos,nbuff;

/*****************************************************************************/
 void par_init_put (char*buff,int*numbuff) 
{
  ibuff=buff;
  nbuff=*numbuff*sizeof(float);
  ipos=0;
}

/*****************************************************************************/
 void par_send (int*mach,int*msgtype)
{
  int ierr=0;
#if defined (RAMS_MPI)
  ierr=MPI_Send(ibuff,ipos,MPI_PACKED,*mach,*msgtype,MPI_COMM_WORLD);
#endif
  /*printf("par_send - %d %d %d %d\n",*mach,*msgtype,ierr,ipos);*/
  if(ierr<0)printf("Error in par_send - %d %d %d\n",*mach,*msgtype,ierr);
}

/*****************************************************************************/
 void par_put_int (int*iwords,int*numwords)
{
  int ierr=0;
#if defined (RAMS_MPI)
  ierr=MPI_Pack(iwords,*numwords,MPI_INT,ibuff,nbuff,&ipos,MPI_COMM_WORLD);
#endif
  if(ierr<0)printf("Error in par_put_int - %d %d\n",ierr,ipos);
}

/*****************************************************************************/
 void par_put_float (float*words,int*numwords)
{
  int ierr=0;
#if defined (RAMS_MPI)
  ierr=MPI_Pack(words,*numwords,MPI_FLOAT,ibuff,nbuff,&ipos,MPI_COMM_WORLD);
#endif
  if(ierr<0)printf("Error in par_put_float - %d\n",ierr);
}

/*****************************************************************************/
 void par_put_char (char*words,int*numbytes)
{
  int ierr=0;
#if defined (RAMS_MPI)
  ierr=MPI_Pack(words,*numbytes,MPI_BYTE,ibuff,nbuff,&ipos,MPI_COMM_WORLD);
#endif
  if(ierr<0)printf("Error in par_put_float - %d\n",ierr);
}

/*****************************************************************************/
 void zero_tag_array ()
{
  int i;
  flag_msgtag = 1;
  for(i=0;i<MAX_MSGTAG;i++)
    {
      mpi_msgtags[i] = 0;
    }
}

/*****************************************************************************/
int par_store_tag (MPI_Request msgtag)
{
   int i;
   for(i=0; i < MAX_MSGTAG; i++) {
      if (mpi_msgtags[i] == 0) {
         mpi_msgtags[i]= msgtag;
         return(i);
      }
   }
  printf("par_store_tag error:msgtag,max_msgtag");
  return (i);
}

/*****************************************************************************/
MPI_Request par_retrieve_tag (int fmsgtag)
{
   MPI_Request i;

   /*For Safety*/
   if (fmsgtag >= MAX_MSGTAG) 
      printf("ERROR:par_retrieve_tag:fmsgtag >= MAX_MSGTAG\n");
    
   i = mpi_msgtags[fmsgtag];

   mpi_msgtags[fmsgtag] = 0;

   return (i);
}

/*****************************************************************************/
 void par_send_noblock (int*mach,int*msgtype,int*fmsgtag)
{
   int ierr=0;
#if defined (RAMS_MPI)
   MPI_Request msgtag;
   /*printf("par_send - %d %d %d %d\n",*mach,*msgtype,ierr,ipos);*/
   if (flag_msgtag == 0) zero_tag_array();
   ierr=MPI_Isend(ibuff,ipos,MPI_PACKED,*mach,*msgtype,MPI_COMM_WORLD,&msgtag);
   *fmsgtag = par_store_tag(msgtag);
   /*printf("par_sent - %d %d %d %d\n",*mach,*msgtype,ierr,ipos);*/
#endif
   /*printf("par_send - %d %d %d %d\n",*mach,*msgtype,ierr,ipos);*/
   if(ierr<0)printf("Error par_send_noblock - %d %d %d \n",*mach,*msgtype,ierr);
}

/*****************************************************************************/
 void par_get_noblock (void*buff,int*numbuff,int*mmtype,int*ihostnum,int*fmsgtag)
{
   int ierr=0;
#if defined (RAMS_MPI)
   MPI_Request msgtag;
   if (flag_msgtag == 0) zero_tag_array();
   nbuff=*numbuff*sizeof(float);
   ipos=0;
   ierr=MPI_Irecv(buff,nbuff,MPI_PACKED,*ihostnum,*mmtype,MPI_COMM_WORLD,&msgtag);
   *fmsgtag = par_store_tag(msgtag);
#endif
   if(ierr<0)printf("Error in par_get_noblock\n");
}

/*****************************************************************************/
 void par_assoc_buff (void*buff,int*numbuff)
{
  ibuff=buff;
  nbuff=*numbuff*sizeof(float);
  ipos=0;
}

/*****************************************************************************/
 void par_wait (int*fmsgtag,int*ibytes,int*msgtype,int*ihostnum)
{
   int ierr=0;
#if defined (RAMS_MPI)
   MPI_Request msgtag;
   MPI_Status status;   
   /*printf("Node waiting - %d %d %d %d\n",*msgtype,*ihostnum,*ibytes,mynum);*/
   if (flag_msgtag == 0) zero_tag_array();
   msgtag = par_retrieve_tag(*fmsgtag);
   ierr=MPI_Wait(&msgtag,&status);
   MPI_Get_count(&status,MPI_PACKED,ibytes);
   *msgtype=status.MPI_TAG;
   *ihostnum=status.MPI_SOURCE;
#endif
   /*printf("Node done wait- %d %d %d %d\n",*msgtype,*ihostnum,*ibytes,mynum);*/
   if(ierr<0)printf("Error in par_wait\n");
}

/*****************************************************************************/
 void par_get_new (void*buff,int*numbuff,int*mmtype,int*ibytes,int*msgtype,int*ihostnum)
{
  int ierr=0;
#if defined (RAMS_MPI)
  MPI_Status status;
  ibuff=buff;
  nbuff=*numbuff*sizeof(float);
  ipos=0;
  /*printf("Node waiting for - %d\n",*mmtype);*/
  ierr=MPI_Recv(ibuff,nbuff,MPI_PACKED,MPI_ANY_SOURCE,*mmtype,MPI_COMM_WORLD,&status);
  MPI_Get_count(&status,MPI_PACKED,ibytes);
  *msgtype=status.MPI_TAG;
  *ihostnum=status.MPI_SOURCE;
#endif
  /*printf("Node got - %d %d %d %d\n",*msgtype,*ihostnum,*ibytes,mynum);*/
  if(ierr<0)printf("Error in par_get_new\n");
}

/*****************************************************************************/
 void par_get_int (int*iwords,int*numwords)
{
  int ierr=0;
#if defined (RAMS_MPI)
  /*printf("Unpack int- %d %d \n",*numwords,nbuff);*/
  ierr=MPI_Unpack(ibuff,nbuff,&ipos,iwords,*numwords,MPI_INT,MPI_COMM_WORLD);
#endif
  /*printf("Node got - %d %d %d %d\n",*msgtype,*ihostnum,*ibytes,mynum);*/
  if(ierr<0)printf("Error in par_get_int-%d \n",ierr);
}

/*****************************************************************************/
 void par_get_float (float*words,int*numwords)
{
  int ierr=0;
#if defined (RAMS_MPI)
  /*printf("Unpack flt- %d %d \n",*numwords,nbuff);*/
  ierr=MPI_Unpack(ibuff,nbuff,&ipos,words,*numwords,MPI_FLOAT,MPI_COMM_WORLD);
#endif
  if(ierr<0)printf("Error in par_get_float-%d \n",ierr);
}

/*****************************************************************************/
 void par_get_char (char*words,int*numbytes)
{
  int ierr=0;
#if defined (RAMS_MPI)
  /*printf("Unpack flt- %d %d \n",*numwords,nbuff);*/
  ierr=MPI_Unpack(ibuff,nbuff,&ipos,words,*numbytes,MPI_BYTE,MPI_COMM_WORLD);
#endif
  if(ierr<0)printf("Error in par_get_float-%d \n",ierr);
}

/*****************************************************************************/
 void par_init_fortran (int*argc,char*fargv,int*farglen,int*machnum,int*machsize)
{
  int i,numarg,carglen;
  char *argvp[20];
  char **argv;

  numarg=*argc;
  carglen=*farglen;
  /*printf("par init numargs: %d %s %d %d\n",numarg,fargv,carglen,*machnum);*/

  for (i = 0; i < numarg; i++) {
    argvp[i]=&(fargv[i*carglen]);
    /*printf("par init args: %d %s %s\n",i,"argvp[i]",argvp[i]);*/
  }
  argv=&(argvp[0]);

#if defined (RAMS_MPI)
  /*printf("par init RAMS_MPI defined \n");*/
  MPI_Init(&numarg, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD,machnum);
  MPI_Comm_size(MPI_COMM_WORLD,machsize);
#endif
  /*printf("par_init: %d %d\n",*machnum,*machsize);*/
}

/*****************************************************************************/
 void par_exit ()
{
#if defined (RAMS_MPI)
  MPI_Finalize();
#endif
/*printf("MP exiting\n");*/
}

/*****************************************************************************/
 void par_pause (int*machnum,int*ibarrier)
{
  int ierr=0;
#if defined (RAMS_MPI)
  ierr=MPI_Barrier(MPI_COMM_WORLD);
#endif
  if(ierr<0)printf("Error in par_pause- %d %d %d\n",*machnum,*ibarrier,ierr);
}

/*****************************************************************************/
 void par_ready (int*nmach,int*machnum,int*ibarrier)
{
  int ierr=0;
#if defined (RAMS_MPI)
  /*printf("par_ready - %d %d\n",*ibarrier,*machnum);*/
  ierr=MPI_Barrier(MPI_COMM_WORLD);
#endif
  if(ierr<0)printf("Error in par_ready - %d %d %d\n",*machnum,*ibarrier,ierr);
}
