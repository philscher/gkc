/*
 * =====================================================================================
 *
 *       Filename:  LennardBernstein.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/18/2012 03:19:23 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */



class CollisionsLennardBernstein : public Collisions

             

	    const cmplxd dfs_dv   = (8.  *(fs[s][m][z][y_k][x][v+1] - fs[s][m][z][y_k][x][v-1]) - (fs[s][m][z][y_k][x][v+2] - fs[s][m][z][y_k][x][v-2]))/(12.*dv);
        const cmplxd ddfs_dvv = (16. *(fs[s][m][z][y_k][x][v+1] + fs[s][m][z][y_k][x][v-1]) - (fs[s][m][z][y_k][x][v+2] + fs[s][m][z][y_k][x][v-2]) - 30.*fs[s][m][z][y_k][x][v])/(12.*pow2(dv));
                                   + collisionBeta  * (g  + V[v] * dfs_dv + v2_rms * ddfs_dvv)
