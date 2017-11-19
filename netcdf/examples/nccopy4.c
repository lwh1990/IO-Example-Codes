#include <stdlib.h>
#include <stdio.h>
#ifndef _WIN32
#include <unistd.h>
#endif
#include <string.h>
#include <netcdf.h>

#define SAME_AS_INPUT (-1)	/* default, if kind not specified */

#define CHECK(stat,f) if(stat != NC_NOERR) {check(stat,#f,__FILE__,__LINE__);} else {}

/* These are in unistd.h; for use with getopt() */
extern int optind;
extern int opterr;
extern char *optarg;

static char *progname; 		/* for error messages */

static void
check(int err, const char* fcn, const char* file, const int line)
{
    fprintf(stderr,"%s\n",nc_strerror(err));
    fprintf(stderr,"Location: function %s; file %s; line %d\n",
	    fcn,file,line);
    fflush(stderr); fflush(stdout);
    exit(1);
}

/* Check error return from malloc, and allow malloc(0) with subsequent free */
static void *
emalloc (size_t size)
{
    void   *p;

    p = (void *) malloc (size==0 ? 1 : size); /* don't malloc(0) */
    if (p == 0) {
	fprintf(stderr,"Out of memory\n");
    }
    return p;
}

/* Forward declaration, because copy_type, copy_vlen_type call each other */
static int copy_type(int igrp, nc_type typeid, int ogrp);

/* get group id in output corresponding to group igrp in input,
 * given parent group id (or root group id) parid in output. */
static int
get_grpid(int igrp, int parid, int *ogrpp) {
    int stat = NC_NOERR;
    int ogid;			/* like igrp but in output file */
    int inparid;

    /* if not root group, get corresponding output groupid from group name */
    stat = nc_inq_grp_parent(igrp, &inparid);
    if(stat == NC_NOERR) {	/* not root group */
	char grpname[NC_MAX_NAME + 1];
	stat = nc_inq_grpname(igrp, grpname);
	CHECK(stat, nc_inq_grpname);
	stat = nc_inq_grp_ncid(parid, grpname, &ogid);
	CHECK(stat, nc_inq_grp_ncid);
    } else if(stat == NC_ENOGRP) { /* root group */
	ogid = parid;
	stat = NC_NOERR;
    } else {
	CHECK(stat, nc_inq_grp_parent);
    }
    *ogrpp = ogid;
    return stat;
}

/* 
 * copy a user-defined variable length type in the group igrp to the
 * group ogrp
 */
static int
copy_vlen_type(int igrp, nc_type itype, int ogrp)
{
    int stat = NC_NOERR; 
    nc_type ibasetype;
    nc_type obasetype;		/* base type in target group */
    char name[NC_MAX_NAME];
    size_t size;
    char basename[NC_MAX_NAME];
    size_t basesize;
    nc_type vlen_type;

    stat = nc_inq_vlen(igrp, itype, name, &size, &ibasetype);
    CHECK(stat, nc_inq_vlen);
    /* to get base type id in target group, use name of base type in
     * source group */
    stat = nc_inq_type(igrp, ibasetype, basename, &basesize);
    CHECK(stat, nc_inq_type);
    stat = nc_inq_typeid(ogrp, basename, &obasetype);
    /* if no such type, create it now */
    if(stat == NC_EBADTYPE) {
	copy_type(igrp, ibasetype, ogrp);
	CHECK(stat, copy_type);
	stat = nc_inq_typeid(ogrp, basename, &obasetype);
    }
    CHECK(stat, nc_inq_typeid);

    /* Now we know base type exists in output and we know its type id */
    stat = nc_def_vlen(ogrp, name, obasetype, &vlen_type);
    CHECK(stat, nc_copy_vlen_type);

    return stat;
}

/* 
 * copy a user-defined opaque type in the group igrp to the group ogrp
 */
static int
copy_opaque_type(int igrp, nc_type itype, int ogrp)
{
    int stat = NC_NOERR; 
    nc_type otype;
    char name[NC_MAX_NAME];
    size_t size;

    stat = nc_inq_opaque(igrp, itype, name, &size);
    CHECK(stat, nc_inq_opaque_type);
    stat = nc_def_opaque(ogrp, size, name, &otype);
    CHECK(stat, copy_opaque_type);

    return stat;
}

/* 
 * copy a user-defined enum type in the group igrp to the group ogrp
 */
static int
copy_enum_type(int igrp, nc_type itype, int ogrp)
{
    int stat = NC_NOERR; 
    nc_type otype;
    nc_type basetype;
    size_t basesize;
    size_t nmembers;
    char name[NC_MAX_NAME];
    int i;

    stat = nc_inq_enum(igrp, itype, name, &basetype, &basesize, &nmembers);
    CHECK(stat, nc_inq_enum);
    stat = nc_def_enum(ogrp, basetype, name, &otype);
    CHECK(stat, nc_def_enum);
    for(i = 0; i < nmembers; i++) { /* insert enum members */
	char ename[NC_MAX_NAME];
	long long val;		/* large enough to hold any integer type */
	stat = nc_inq_enum_member(igrp, itype, i, ename, &val);
	CHECK(stat, nc_inq_enum_member);
	stat = nc_insert_enum(ogrp, otype, ename, &val);
	CHECK(stat, nc_insert_enum);
    }
    return stat;
}

/* 
 * copy a user-defined compound type in the group igrp to the group ogrp
 */
static int
copy_compound_type(int igrp, nc_type itype, int ogrp)
{
    int stat = NC_NOERR; 
    char name[NC_MAX_NAME];
    size_t size;
    size_t nfields;
    nc_type otype;
    int fid;

    stat = nc_inq_compound(igrp, itype, name, &size, &nfields);
    CHECK(stat, nc_inq_compound);
    stat = nc_def_compound(ogrp, size, name, &otype);
    CHECK(stat, nc_def_compound);

    for (fid = 0; fid < nfields; fid++) {
	char fname[NC_MAX_NAME];
	char ftypename[NC_MAX_NAME];
	size_t foff;
	nc_type iftype, oftype;
	int fndims, fdimsizes[NC_MAX_DIMS];

	stat = nc_inq_compound_field(igrp, itype, fid, fname, &foff, &iftype, 
				     &fndims, fdimsizes);
	CHECK(stat, nc_inq_compound_field);
	/* type ids in source don't necessarily correspond to same
	 * typeids in destination, so look up destination typeid by using
	 * field type name */
	stat = nc_inq_type(igrp, iftype, ftypename, NULL);
	CHECK(stat, nc_inq_type);
	stat = nc_inq_typeid(ogrp, ftypename, &oftype);
	CHECK(stat, nc_inq_typeid);
	if(fndims == 0) {
	    stat = nc_insert_compound(ogrp, otype, fname, foff, oftype);
	    CHECK(stat, nc_insert_compound);	    
	} else {		/* field is array type */
	    stat = nc_insert_array_compound(ogrp, otype, fname, foff, 
					    oftype, fndims, fdimsizes);
	    CHECK(stat, nc_insert_array_compound);
	}
    }
    return stat;
}


/* 
 * copy a user-defined type in the group igrp to the group ogrp
 */
static int
copy_type(int igrp, nc_type typeid, int ogrp)
{
    int stat = NC_NOERR; 
    nc_type type_class;

    stat = nc_inq_user_type(igrp, typeid, NULL, NULL, NULL, NULL, &type_class); 
    CHECK(stat, nc_inq_user_type);

    switch(type_class) {
    case NC_VLEN:
	stat = copy_vlen_type(igrp, typeid, ogrp);
	CHECK(stat, copy_vlen_type);
	break;
    case NC_OPAQUE:
	stat = copy_opaque_type(igrp, typeid, ogrp);
	CHECK(stat, copy_opaque_type);
	break;
    case NC_ENUM:
	stat = copy_enum_type(igrp, typeid, ogrp);
	CHECK(stat, copy_enum_type);
	break;
    case NC_COMPOUND:
	stat = copy_compound_type(igrp, typeid, ogrp);
	CHECK(stat, copy_compound_type);
	break;
    default:
	stat = NC_EBADTYPE;
	CHECK(stat, copy_type);
    }
    return stat;
}

/* Copy a group and all its subgroups, recursively, from group igrp in
 * input to parent group ogrp in destination.  This just creates all
 * the groups in the destination, but doesn't copy anything that's in
 * the groups. */
static int
copy_groups(int igrp, int ogrp)
{
    int stat = NC_NOERR;
    int inparid;
    int ogid;			/* like igrp but in output file */
    int numgrps;
    int *grpids;
    int i;

    /* if not root group, create corresponding new group in ogrp */
    stat = nc_inq_grp_parent(igrp, &inparid);
    if(stat == NC_NOERR) {
	/* create new subgroup */
	char grpname[NC_MAX_NAME + 1];
	stat = nc_inq_grpname(igrp, grpname);
	CHECK(stat, nc_inq_grpname);
	stat = nc_def_grp(ogrp, grpname, &ogid);
	CHECK(stat, nc_def_grp);
    } else if(stat == NC_ENOGRP) {
	ogid = ogrp;
	stat = NC_NOERR;
    } else {
	CHECK(stat, nc_inq_grp_parent);
    }
    
    /* Copy any subgroups */
    stat = nc_inq_grps(igrp, &numgrps, NULL);
    grpids = (int *)emalloc(sizeof(int) * (numgrps + 1));
    stat = nc_inq_grps(igrp, &numgrps, grpids);
    CHECK(stat, nc_inq_grps);

    for(i = 0; i < numgrps; i++) {
	stat = copy_groups(grpids[i], ogid);
	CHECK(stat, copy_group);
    }
    free(grpids);
    return stat;    
}

/* 
 * Copy the user-defined types in this group (igrp) and all its
 * subgroups, recursively, to corresponding group in output (ogrp)
 */
static int
copy_types(int igrp, int ogrp)
{
    int stat = NC_NOERR; 
    int ntypes;
    nc_type *types = NULL;
    int numgrps;
    int *grpids = NULL;
    int i;

    stat = nc_inq_typeids(igrp, &ntypes, NULL);
    CHECK(stat, nc_inq_typeids);

    if(ntypes > 0) {
	types = (nc_type *) emalloc(ntypes * sizeof(nc_type));
	stat = nc_inq_typeids(igrp, &ntypes, types);
	CHECK(stat, nc_inq_typeids);
	for (i = 0; i < ntypes; i++) {
	    stat = copy_type(igrp, types[i], ogrp);
	    CHECK(stat, copy_type);
	}
	free(types);
    }

    /* Copy types from subgroups */
    stat = nc_inq_grps(igrp, &numgrps, NULL);
    CHECK(stat, nc_inq_grps);
    if(numgrps > 0) {
	grpids = (int *)emalloc(sizeof(int) * (numgrps + 1));
	stat = nc_inq_grps(igrp, &numgrps, grpids);
	CHECK(stat, nc_inq_grps);
	for(i = 0; i < numgrps; i++) {
	    int ogid;
	    /* get groupid in output corresponding to grpids[i] in
	     * input, given parent group (or root group) ogrp in
	     * output */
	    stat = get_grpid(grpids[i], ogrp, &ogid);
	    CHECK(stat, get_grpid);
	    stat = copy_types(grpids[i], ogid);
	    CHECK(stat, copy_types);
	}
	free(grpids);
    }
    return stat;
}

/* Copy all netCDF-4 specific variable properties such as chunking,
 * endianness, deflation, checksumming, fill, etc. */
static int
copy_var_specials(int igrp, int varid, int ogrp, int o_varid)
{
    int stat = NC_NOERR;
    {				/* handle chunking parameters */
	int ndims;
	stat = nc_inq_varndims(igrp, varid, &ndims);
	CHECK(stat, nc_inq_varndims);
	if (ndims > 0) {		/* no chunking for scalar variables */
	    int contig = 0;
	    stat = nc_inq_var_chunking(igrp, varid, &contig, NULL);
	    CHECK(stat, nc_inq_var_chunking);
	    if(contig == 1) {
		stat = nc_def_var_chunking(ogrp, o_varid, NC_CONTIGUOUS, NULL);
		CHECK(stat, nc_def_var_chunking);
	    } else {
		size_t chunkp[NC_MAX_DIMS];
		stat = nc_inq_var_chunking(igrp, varid, NULL, chunkp);
		CHECK(stat, nc_inq_var_chunking);
		/* explicitly set chunking, even if default */
		stat = nc_def_var_chunking(ogrp, o_varid, NC_CHUNKED, chunkp);
		CHECK(stat, nc_def_var_chunking);
	    }
	}
    }
    {				/* handle compression parameters */
	int shuffle=NC_NOSHUFFLE, deflate=0, deflate_level=0;
	stat = nc_inq_var_deflate(igrp, varid, 
				  &shuffle, &deflate, &deflate_level);
	CHECK(stat, nc_inq_var_deflate);
	if(deflate != 0 || shuffle != 0) {
	    stat = nc_def_var_deflate(ogrp, o_varid, 
				      shuffle, deflate, deflate_level);
	    CHECK(stat, nc_def_var_deflate);
	}
    }
    {				/* handle checksum parameters */
	int fletcher32 = 0;
	stat = nc_inq_var_fletcher32(igrp, varid, &fletcher32);
	CHECK(stat, nc_inq_var_fletcher32);
	if(fletcher32 != 0) {
	    stat = nc_def_var_fletcher32(ogrp, o_varid, fletcher32);
	    CHECK(stat, nc_def_var_fletcher32);
	}
    }
    {				/* handle endianness */
	int endianness = 0;
	stat = nc_inq_var_endian(igrp, varid, &endianness);
	CHECK(stat, nc_inq_var_endian);
	if(endianness != NC_ENDIAN_NATIVE) { /* native is the default */
	    stat = nc_def_var_endian(ogrp, o_varid, endianness);
	    CHECK(stat, nc_def_var_endian);
	}
    }
    return stat;
}

/* Free an array of vlens.  This might be useful in the library. */
static int
nc_free_vlens(
    size_t len,			/* number of vlens in array to free */
    nc_vlen_t vlens[]		/* array of vlens */
    ) {
    int stat = NC_NOERR;
    size_t i;
    for(i = 0; i < len; i++) {
    	stat = nc_free_vlen(&vlens[i]);
    	CHECK(stat, nc_free_vlen);
    }
    return stat;
}

/* Copy dimensions from group igrp to group ogrp */
static int
copy_dims(int igrp, int ogrp)
{
    int stat = NC_NOERR;
    int ndims;
    int nunlims;
    int dgrp;
    int dimids[NC_MAX_DIMS];
    int unlimids[NC_MAX_DIMS];

    stat = nc_inq_ndims(igrp, &ndims);
    CHECK(stat, nc_inq_ndims);

   /* In netCDF-4 files, dimids may not be sequential because they
    * may be defined in various groups, and we are only looking at one
    * group at a time. */
    /* Find the dimension ids in this group, don't include parents. */
    stat = nc_inq_dimids(igrp, NULL, dimids, 0);
    CHECK(stat, nc_inq_dimids);
    /* Find the number of unlimited dimensions and get their IDs */
    stat = nc_inq_unlimdims(igrp, &nunlims, unlimids);
    CHECK(stat, nc_inq_unlimdims);

    /* Copy each dimension to output, including unlimited dimension(s) */
    for (dgrp = 0; dgrp < ndims; dgrp++) {
	char name[NC_MAX_NAME];
	size_t length;
	int is_unlim;
	int uld;
	int dimid;

	is_unlim = 0;
	dimid = dimids[dgrp];
	for (uld = 0; uld < nunlims; uld++) {
	    if(dimid == unlimids[uld]) {
		is_unlim = 1;
		break;
	    }	  
	}

	stat = nc_inq_dim(igrp, dimid, name, &length);
	if (stat == NC_EDIMSIZE && sizeof(size_t) < 8) {
	    fprintf(stderr, "dimension \"%s\" requires 64-bit platform\n", 
		    name);
	}	
	CHECK(stat, nc_inq_dim);
	if(is_unlim) {
	    stat = nc_def_dim(ogrp, name, NC_UNLIMITED, NULL);
	} else {
	    stat = nc_def_dim(ogrp, name, length, NULL);
	}
	CHECK(stat, nc_def_dim);	
    }
    return stat;
}

/* Copy the attributes for variable ivar in group igrp to variable
 * ovar in group ogrp.  Global (group) attributes are specified by
 * using the varid NC_GLOBAL */
static int
copy_atts(int igrp, int ivar, int ogrp, int ovar)
{
    int natts;
    int iatt;
    int stat = NC_NOERR;

    stat = nc_inq_varnatts(igrp, ivar, &natts);
    CHECK(stat, nc_inq_varnatts);
    
    for(iatt = 0; iatt < natts; iatt++) {
	char name[NC_MAX_NAME];
	stat = nc_inq_attname(igrp, ivar, iatt, name);
	CHECK(stat, nc_inq_attname);
	stat = nc_copy_att(igrp, ivar, name, ogrp, ovar);
	CHECK(stat, nc_copy_att);
    }
    return stat;
}

/* copy the schema for a single variable in group igrp (from a file of
 * kind inkind) to group ogrp */
static int
copy_var(int igrp, int inkind, int varid, int ogrp)
{
    int stat = NC_NOERR;
    int ndims;
    int idimids[NC_MAX_DIMS];	/* ids of dims for input variable */
    int odimids[NC_MAX_DIMS];	/* ids of dims for output variable */
    char name[NC_MAX_NAME];
    nc_type typeid, o_typeid;
    int natts;
    int i;
    int o_varid;

    stat = nc_inq_varndims(igrp, varid, &ndims);
    CHECK(stat, nc_inq_varndims);
    stat = nc_inq_var(igrp, varid, name, &typeid, &ndims, idimids, &natts);
    CHECK(stat, nc_inq_var);
    o_typeid = typeid;
    if (typeid > NC_STRING) {	/* user-defined type */
	/* type ids in source don't necessarily correspond to same
	 * typeids in destination, so look up destination typeid by
	 * using type name */
	char type_name[NC_MAX_NAME];
	stat = nc_inq_type(igrp, typeid, type_name, NULL);
	CHECK(stat, nc_inq_type);
	stat = nc_inq_typeid(ogrp, type_name, &o_typeid);
	CHECK(stat, nc_inq_typeid);
    }

    /* get the corresponding dimids in the output file */
    for(i = 0; i < ndims; i++) {
	char dimname[NC_MAX_NAME];
	stat = nc_inq_dimname(igrp, idimids[i], dimname);
	CHECK(stat, nc_inq_dimname);
	stat = nc_inq_dimid(ogrp, dimname, &odimids[i]);
	CHECK(stat, nc_inq_dimid);
    }

    /* define the output variable */
    stat = nc_def_var(ogrp, name, o_typeid, ndims, odimids, &o_varid);
    CHECK(stat, nc_def_var);

    /* attach the variable attributes to the output variable */
    stat = copy_atts(igrp, varid, ogrp, o_varid);
    CHECK(stat, copy_atts);
    if(inkind == NC_FORMAT_NETCDF4 || inkind == NC_FORMAT_NETCDF4_CLASSIC) {
	/* Copy all netCDF-4 specific variable properties such as
	 * chunking, endianness, deflation, checksumming, fill, etc. */
	stat = copy_var_specials(igrp, varid, ogrp, o_varid);
	CHECK(stat, copy_var_specials);
    }
    return stat;
}

/* copy the schema for all the variables in group igrp (from a file of kind
 * inkind) to group ogrp */
static int
copy_vars(int igrp, int inkind, int ogrp)
{
    int stat = NC_NOERR;
    int nvars;
    int varid;
    
    stat = nc_inq_nvars(igrp, &nvars);
    CHECK(stat, nc_inq_nvars);
    for (varid = 0; varid < nvars; varid++) {
	stat = copy_var(igrp, inkind, varid, ogrp);
	CHECK(stat, copy_var);
	
    }
    return stat;
}

/* Copy the schema in a group and all its subgroups, recursively, from
 * group igrp in input to parent group ogrp in destination.  inkind
 * is the kind of netCDF file that igrp identifies (1 -> classic, 2 -<
 * 64-bit offset, 3-> netCDF-4, 4-> netCDF-4 classic model). */
static int
copy_schema(int igrp, int inkind, int ogrp) 
{
    int stat = NC_NOERR;
    int ogid;			/* like igrp but in output file */
    int i;

    /* get groupid in output corresponding to group igrp in input,
     * given parent group (or root group) ogrp in output */
    stat = get_grpid(igrp, ogrp, &ogid);
    CHECK(stat, get_grpid);

    stat = copy_dims(igrp, ogid);
    CHECK(stat, copy_dims);
    stat = copy_atts(igrp, NC_GLOBAL, ogid, NC_GLOBAL);
    CHECK(stat, copy_atts);
    stat = copy_vars(igrp, inkind, ogid);
    CHECK(stat, copy_vars);
    {
	int numgrps;
	int *grpids;
	/* Copy schema from subgroups */
	stat = nc_inq_grps(igrp, &numgrps, NULL);
	grpids = (int *)emalloc(sizeof(int) * (numgrps + 1));
	stat = nc_inq_grps(igrp, &numgrps, grpids);
	CHECK(stat, nc_inq_grps);
	
	for(i = 0; i < numgrps; i++) {
	    stat = copy_schema(grpids[i], inkind, ogid);
	    CHECK(stat, copy_schema);
	}
	free(grpids);
    }
    return stat;    
}

/* Return number of values for a variable varid in a group igrp, as
 * well as start and count arrays for array access, assumed to be
 * preallocated to hold one value for each dimension of variable */
static int
inq_nvals(int igrp, int varid, 
	  size_t *startp, size_t *countp, size_t *nvalsp) {
    int stat = NC_NOERR;
    int ndims;
    int dimids[NC_MAX_DIMS];
    int dim;
    size_t nvals = 1;

    stat = nc_inq_varndims(igrp, varid, &ndims);
    CHECK(stat, nc_inq_varndims);
    stat = nc_inq_vardimid (igrp, varid, dimids);
    CHECK(stat, nc_inq_vardimid);
    for(dim = 0; dim < ndims; dim++) {
	size_t len;
	stat = nc_inq_dimlen(igrp, dimids[dim], &len);
	CHECK(stat, nc_inq_dimlen);
	nvals *= len;
	startp[dim] = 0;
	countp[dim] = len;
    }
    *nvalsp = nvals;
    return stat;
}

/* From netCDF type in group igrp, get size in memory needed for each
 * value */
static int
inq_value_size(int igrp, nc_type vartype, size_t *value_sizep) {
    int stat = NC_NOERR;
    
    stat = nc_inq_type(igrp, vartype, NULL, value_sizep);
    CHECK(stat, nc_inq_type);
    return stat;
}

/* Copy data from variable varid in group igrp to corresponding group
 * ogrp */
static int
copy_var_data(int igrp, int inkind, int varid, int ogrp) {
    int stat = NC_NOERR;
    nc_type vartype;
    size_t value_size;		/* size in bytes of one value of variable */
    size_t nvalues;		/* number of values for this variable */
    void *buf;			/* buffer for the variable values */
    char varname[NC_MAX_NAME];
    int ovarid;
    size_t start[NC_MAX_DIMS];
    size_t count[NC_MAX_DIMS];

    /* Note: for simplicity, this example code copies a whole variable
     * at a time, so all the data's variable must fit in memory.  This
     * won't work for variables too large to fit in memory.  For an
     * approach that works efficiently in the general case for any
     * size variable, see the use of nc_get_iter() and nc_next_iter()
     * in the nccopy utlity that comes with the netCDF
     * distribution.  */
    stat = inq_nvals(igrp, varid, start, count, &nvalues);
    CHECK(stat, inq_nvals);
    if(nvalues == 0)
	return stat;
    stat = nc_inq_vartype(igrp, varid, &vartype);
    CHECK(stat, nc_inq_vartype);
    /* from type, get size in memory needed for each value */
    stat = inq_value_size(igrp, vartype, &value_size);
    CHECK(stat, inq_value_size);
    /* get corresponding output variable */
    stat = nc_inq_varname(igrp, varid, varname);
    CHECK(stat, nc_inq_varname);
    stat = nc_inq_varid(ogrp, varname, &ovarid);
    CHECK(stat, nc_inq_varid);
    buf = emalloc(value_size * nvalues);
    if(inkind ==  NC_FORMAT_NETCDF4) {
	stat = nc_get_vara(igrp, varid, start, count, buf);
	CHECK(stat, nc_get_vara);
	stat = nc_put_vara(ogrp, ovarid, start, count, buf);
	CHECK(stat, nc_put_vara);
	/* we have to explicitly free values for strings and vlens */
	if(vartype == NC_STRING) {
	    stat = nc_free_string(nvalues, (char **)buf);
	    CHECK(stat, nc_free_string);
	} else if(vartype > NC_STRING) { /* user-defined type */
	    nc_type vclass;
	    stat = nc_inq_user_type(igrp, vartype, NULL, NULL, NULL, 
				    NULL, &vclass); 
	    CHECK(stat, nc_inq_user_type);
	    if(vclass == NC_VLEN) {
		stat = nc_free_vlens(nvalues, (nc_vlen_t *)buf);
		CHECK(stat, nc_free_vlens);
	    }
	}
    } else {
	/* Unfortunately, above typeless copy not allowed for
	 * classic model */
	switch(vartype) {
	case NC_BYTE:
	    stat = nc_get_vara_schar(igrp, varid, start, count, buf);
	    CHECK(stat, nc_get_vara_schar);
	    stat = nc_put_vara_schar(ogrp, ovarid, start, count, buf);
	    CHECK(stat, nc_put_vara_schar);
	    break;
	case NC_CHAR:
	    stat = nc_get_vara_text(igrp, varid, start, count, buf);
	    CHECK(stat, nc_get_vara_text);
	    stat = nc_put_vara_text(ogrp, ovarid, start, count, buf);
	    CHECK(stat, nc_put_vara_text);
	    break;
	case NC_SHORT:
	    stat = nc_get_vara_short(igrp, varid, start, count, buf);
	    CHECK(stat, nc_get_vara_short);
	    stat = nc_put_vara_short(ogrp, ovarid, start, count, buf);
	    CHECK(stat, nc_put_vara_short);
	    break;
	case NC_INT:
	    stat = nc_get_vara_int(igrp, varid, start, count, buf);
	    CHECK(stat, nc_get_vara_int);
	    stat = nc_put_vara_int(ogrp, ovarid, start, count, buf);
	    CHECK(stat, nc_put_vara_int);
	    break;
	case NC_FLOAT:
	    stat = nc_get_vara_float(igrp, varid, start, count, buf);
	    CHECK(stat, nc_get_vara_float);
	    stat = nc_put_vara_float(ogrp, ovarid, start, count, buf);
	    CHECK(stat, nc_put_vara_float);
	    break;
	case NC_DOUBLE:
	    stat = nc_get_vara_double(igrp, varid, start, count, buf);
	    CHECK(stat, nc_get_vara_double);
	    stat = nc_put_vara_double(ogrp, ovarid, start, count, buf);
	    CHECK(stat, nc_put_vara_double);
	    break;
	default:
	    CHECK(NC_EBADTYPE, copy_var_data);
	    break;
	}
    }
    free(buf);
    return stat;
}

/* Copy data from variables in group igrp to variables in
 * corresponding group with parent ogrp, and all subgroups
 * recursively  */
static int
copy_data(int igrp, int inkind, int ogrp)
{
    int stat = NC_NOERR;
    int ogid;
    int numgrps;
    int *grpids;
    int i;
    int nvars;
    int varid;

    /* get groupid in output corresponding to group igrp in input,
     * given parent group (or root group) ogrp in output */
    stat = get_grpid(igrp, ogrp, &ogid);
    CHECK(stat, get_grpid);
    
    /* Copy data from this group */
    stat = nc_inq_nvars(igrp, &nvars);
    CHECK(stat, nc_inq_nvars);
    for (varid = 0; varid < nvars; varid++) {
	stat = copy_var_data(igrp, inkind, varid, ogid);
	CHECK(stat, copy_var_data);
    }
    /* Copy data from subgroups */
    stat = nc_inq_grps(igrp, &numgrps, NULL);
    grpids = (int *)emalloc(sizeof(int) * (numgrps + 1));
    stat = nc_inq_grps(igrp, &numgrps, grpids);
    CHECK(stat, nc_inq_grps);

    for(i = 0; i < numgrps; i++) {
	stat = copy_data(grpids[i], inkind, ogid);
	CHECK(stat, copy_data);
    }
    free(grpids);
    return stat;
}

/* copy infile to outfile using netCDF API, kind specifies which
 * netCDF format for output: 0 -> same as input, 1 -> classic, 2 ->
 * 64-bit offset, 3 -> netCDF-4, 4 -> netCDF-4 classic model */
static int
copy(char* infile, char* outfile, int kind)
{
    int stat = NC_NOERR;
    int igrp, ogrp;
    int inkind, outkind;

    stat = nc_open(infile,NC_NOWRITE,&igrp);
    CHECK(stat,nc_open);

    stat = nc_inq_format(igrp, &inkind);
    CHECK(stat,nc_inq_format);

    outkind = (kind == SAME_AS_INPUT) ? inkind : kind;

    switch(outkind) {
    case NC_FORMAT_CLASSIC:
	stat = nc_create(outfile,NC_CLOBBER,&ogrp);
	break;
    case NC_FORMAT_64BIT:
	stat = nc_create(outfile,NC_CLOBBER|NC_64BIT_OFFSET,&ogrp);
	break;
    case NC_FORMAT_NETCDF4:
	stat = nc_create(outfile,NC_CLOBBER|NC_NETCDF4,&ogrp);
	break;
    case NC_FORMAT_NETCDF4_CLASSIC:
	stat = nc_create(outfile,NC_CLOBBER|NC_NETCDF4|NC_CLASSIC_MODEL,&ogrp);
	break;
    default:
	fprintf(stderr,"%s: bad value (%d) for -k option\n", progname, kind);
	exit(1);
    }
    CHECK(stat,nc_create);

    /* Because types in one group may depend on types in a different
     * group, need to create all groups before defining types */
    if(inkind == NC_FORMAT_NETCDF4) { /* enhanced data model with
				       * groups, user-defined types */
	stat = copy_groups(igrp, ogrp);
	CHECK(stat,copy_groups);
	stat = copy_types(igrp, ogrp);
	CHECK(stat,copy_types);
    }
    stat = copy_schema(igrp, inkind, ogrp);
    CHECK(stat,copy_schema);
    stat = nc_enddef(ogrp);
    CHECK(stat, nc_enddef);
    stat = copy_data(igrp, inkind, ogrp);
    CHECK(stat,copy_data);

    stat = nc_close(igrp);
    CHECK(stat,nc_close);
    stat = nc_close(ogrp);
    CHECK(stat,nc_close);
    return stat;
}

static void
usage(void)
{
#define USAGE   "\
  [-k n]    kind of netCDF format for output file, default same as input\n\
	    1 classic, 2 64-bit offset, 3 netCDF-4, 4 netCDF-4 classic model\n\
  infile    name of netCDF input file\n\
  outfile   name for netCDF output file\n"

    (void) fprintf(stderr,
		   "%s [-k n] infile outfile\n%s",
		   progname,
		   USAGE);
}

int
main(int argc, char**argv)
{
    char* inputfile = NULL;
    char* outputfile = NULL;
    int kind = SAME_AS_INPUT; /* default, output same format as input */
    int c;

    opterr = 1;
    progname = argv[0];

    if (argc <= 1)
    {
       usage();
       return 1;
    }

    while ((c = getopt(argc, argv, "k:")) != EOF) {
	switch(c) {
	case 'k':
	    kind = strtol(optarg, NULL, 10);
	    break;
	default: 
	    usage();
	    exit(1);
	    break;
        }
    }

    argc -= optind;
    argv += optind;

    if (argc != 2) {
	fprintf(stderr,"one input file and one output file required\n");
	exit(1);
    }
    inputfile = argv[0];
    outputfile = argv[1];

    if(strcmp(inputfile, outputfile) == 0) {
	fprintf(stderr,"output would overwrite input\n");
	exit(1);
    }

    if(copy(inputfile, outputfile, kind) != NC_NOERR)
        exit(1);
    return 0;
}
