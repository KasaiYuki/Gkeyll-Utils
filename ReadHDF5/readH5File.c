// https://support.hdfgroup.org/ftp/HDF5/examples/examples-by-api/hdf5-examples/1_8/C/H5G/h5ex_g_traverse.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hdf5.h"
#define FILE_NAME "squaregrid.hdf"

typedef struct xy_t {
   float r1;
   float r2;
} xy_t;

char *names[5];
char *groupName;
int iter = 0;

/*
 * Define operator data structure type for H5Literate callback.
 * During recursive iteration, these structures will form a
 * linked list that can be searched for duplicate groups,
 * preventing infinite recursion.
 */
struct opdata {
   unsigned recurs;     /* Recursion level.  0=root */
   struct opdata *prev; /* Pointer to previous opdata */
   haddr_t addr;        /* Group address */
};

/*
 * Operator function to be called by H5Literate.
 */
herr_t op_func(hid_t loc_id, const char *name, const H5L_info_t *info,
               void *operator_data);

/*
 * Function to check for duplicate groups in a path.
 */
int group_check(struct opdata *od, haddr_t target_addr);

int populateStructures(int numDim, int*** domBox, int*** interBox, int*** mask, float**** xy);

int main(void)
{
   setvbuf(stdin, NULL, _IONBF, 0);
   // VALUES TO OUTPUT: ASSUMING 2 PHYSICAL DIMENSIONS OF DATA
   int number_of_dimensions = -1;
   int** domain_box;
   int** interior_box;
   int** mask;
   float*** xy;
   
   int rc = populateStructures(number_of_dimensions, &domain_box, &interior_box, &mask, &xy);
   #ifdef DEBUG
   for(int i = 0; i < 5; i++) {
      for(int j = 0; j < 5; j++) 
         printf("(%f,%f) ", xy[i][j][0], xy[i][j][1]);
      printf("\n");
   }
   #endif
   return 0;
}

/************************************************************

  Operator function.  This function prints the name and type
  of the object passed to it.  If the object is a group, it
  is first checked against other groups in its path using
  the group_check function, then if it is not a duplicate,
  H5Literate is called for that group.  This guarantees that
  the program will not enter infinite recursion due to a
  circular path in the file.

 ************************************************************/
herr_t op_func(hid_t loc_id, const char *name, const H5L_info_t *info,
               void *operator_data)
{
   herr_t status, return_val = 0;
   H5O_info_t infobuf;
   struct opdata *od = (struct opdata *)operator_data;
   /* Type conversion */
   unsigned spaces = 2 * (od->recurs + 1);
   /* Number of whitespaces to prepend
      to output */

   /*
    * Get type of the object and display its name and type.
    * The name of the object is passed to this function by
    * the Library.
    */
   status = H5Oget_info_by_name(loc_id, name, &infobuf, H5P_DEFAULT);
   printf("%*s", spaces, ""); /* Format output */
   switch (infobuf.type) {
      case H5O_TYPE_GROUP:
         printf("Group: %s {\n", name);
         groupName = strdup(name);
         /*
          * Check group address against linked list of operator
          * data structures.  We will always run the check, as the
          * reference count cannot be relied upon if there are
          * symbolic links, and H5Oget_info_by_name always follows
          * symbolic links.  Alternatively we could use H5Lget_info
          * and never recurse on groups discovered by symbolic
          * links, however it could still fail if an object's
          * reference count was manually manipulated with
          * H5Odecr_refcount.
          */
         if (group_check(od, infobuf.addr)) {
            printf("%*s  Warning: Loop detected!\n", spaces, "");
         }
         else {
            /*
             * Initialize new operator data structure and
             * begin recursive iteration on the discovered
             * group.  The new opdata structure is given a
             * pointer to the current one.
             */
            struct opdata nextod;
            nextod.recurs = od->recurs + 1;
            nextod.prev = od;
            nextod.addr = infobuf.addr;
            return_val = H5Literate_by_name(loc_id, name, H5_INDEX_NAME,
                                            H5_ITER_NATIVE, NULL, op_func, (void *)&nextod,
                                            H5P_DEFAULT);
         }
         printf("%*s}\n", spaces, "");
         break;
      case H5O_TYPE_DATASET:
         printf("Dataset: %s\n", name);
         strcpy(names[iter], name);
         break;
      case H5O_TYPE_NAMED_DATATYPE:
         printf("Datatype: %s\n", name);
         break;
      default:
         printf("Unknown: %s\n", name);
   }
   iter++;
   return return_val;
}

/************************************************************

  This function recursively searches the linked list of
  opdata structures for one whose address matches
  target_addr.  Returns 1 if a match is found, and 0
  otherwise.

 ************************************************************/
int group_check(struct opdata *od, haddr_t target_addr)
{
   if (od->addr == target_addr)
      return 1; /* Addresses match */
   else if (!od->recurs)
      return 0; /* Root group reached with no matches */
   else
      return group_check(od->prev, target_addr);
   /* Recursively examine the next node */
}

int populateStructures(int numDim, int ***domBox, int ***interBox, int ***mask, float ****xy)
{
   // ----------------
   hid_t file; 
   herr_t status;
   H5O_info_t infobuf;
   struct opdata od;

   /*
    * Open file and initialize the operator data structure.
    */
   file = H5Fopen(FILE_NAME, H5F_ACC_RDONLY, H5P_DEFAULT);
   status = H5Oget_info(file, &infobuf);
   od.recurs = 0;
   od.prev = NULL;
   od.addr = infobuf.addr;

   printf("/ {\n");
   for (int j = 0; j < 5; ++j)
      names[j] = calloc(1, 100);
   status = H5Literate(file, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, op_func,
                       (void *)&od); // Iterate through datasets in group
   printf("}\n");
   hid_t grp = H5Gopen1(file, groupName);
   for(int i = 0; i < 5; i++) {
      hid_t dset = H5Dopen1(grp, names[i]);
      hid_t dspace = H5Dget_space(dset);
      const int ndims = H5Sget_simple_extent_ndims(dspace); // number of dimensions
      hsize_t *dims = calloc(ndims, sizeof(hsize_t)); // will hold the size of the dataset for each dimension
                                                      // dims[0] = y, dims[1] = x
      status = H5Sget_simple_extent_dims(dspace, dims, NULL);

      #ifdef DEBUG
      for(int i = 0; i < ndims; i++) {
         printf(" %lld |", dims[i]);
      }
      #endif
      if(strcmp(names[i], "dim") == 0) { // ndims = 1
         int (*vals)[dims[0]] = calloc(dims[1], sizeof *vals);
         status = H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vals);
         numDim = vals[0][0];
      }
      else if(strcmp(names[i], "xy") == 0) { // ndims = 3
         float (*vals)[dims[0]][dims[2]] = calloc(dims[1], sizeof *vals);
         // float vals[dims[0]][dims[1]][dims[2]]; // Fully static array version

         status = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vals);
         *xy = calloc(dims[0], sizeof(float**));

         #ifdef DEBUG
         printf("\n");
         for(int j = 0; j < 10; j++) {
            for(int k = 0; k < 10; k++) {
               printf("(%f,%f)\n ", vals[j][k][0], vals[j][k][1]);
            }
            printf("\n");
         }
         #endif
         for(int j = 0; j < dims[0]; j++) {
            (*xy)[j] = calloc(dims[1], sizeof(float*));
            for(int k = 0; k < dims[1]; k++) {
               (*xy)[j][k] = calloc(dims[2], sizeof(float));
               for(int l = 0; l < dims[2]; l++) {
                  (*xy)[j][k][l] = vals[j][k][l];
               }
            }
         }        
      }
      else {
         int (*vals)[dims[0]] = calloc(dims[1], sizeof *vals);
         

         status = H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vals);

         #ifdef DEBUG
         printf("---------%d-------------\n", strcmp(names[i], "domain_box"));
         for(int j = 0; j < dims[0]; j++) {
            for(int k = 0; k < dims[0]; k++) {
               printf("%d ", vals[j][k]);
            }
            printf("\n");
         }
         printf("----------------------\n");
         #endif
         if(strcmp(names[i], "domain_box") == 0) {

            *domBox = calloc(dims[0], sizeof(int*));
            for(int j = 0; j < dims[0]; j++) {
               (*domBox)[j] = calloc(dims[1], sizeof(int));
               for(int k = 0; k < dims[1]; k++) {
                  (*domBox)[j][k] = vals[j][k];
               }
            }
         }
         else if(strcmp(names[i], "interior_box") == 0) {
            *interBox = calloc(dims[0], sizeof(int*));
            for(int j = 0; j < dims[0]; j++) {
               (*interBox)[j] = calloc(dims[1], sizeof(int));
               for(int k = 0; k < dims[1]; k++)
                  (*interBox)[j][k] = vals[j][k];
            }
         }
         else { // mask
            *mask = calloc(dims[0], sizeof(int*));
            for(int j = 0; j < dims[0]; j++) {
               (*mask)[j] = calloc(dims[1], sizeof(int));
               for(int k = 0; k < dims[1]; k++)
                  (*mask)[j][k] = vals[j][k];
            }
         }
         
      }
      free(dims);
   }
   for(int i = 0; i < 5; i++)
      free(names[i]);

   status = H5Fclose(file);
   return 0;
}
