/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-06-24 09:39:53
 * @LastEditTime: 2024-08-10 14:28:13
 * @FilePath: /cpgrid/src/petscdmcp.C
 * @Description: reference from CP
 *
 */
#include "petscdmcp.h"
#include "unstruct_mesh.h"
#include <map>
#include <set>
#include <vector>
// PETSc includes
#define PETSCDM_DLL
#include <petsc/private/dmpleximpl.h> /*I   "petscdmplex.h"   I*/
#include <petsc/private/hashseti.h>   /*I   "petscdmplex.h"   I*/
#include <petscsf.h>
#include <petscdmplextransform.h>
#include <petscdmlabelephemeral.h>
#include <petsc/private/kernels/blockmatmult.h>
#include <petsc/private/kernels/blockinvert.h>
#include <petscdm.h>
// my includes
#include "unstruct_mesh.h"
#include "system.h"
EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "DMCreate_CP"
#define DMCP_NO_DECOMPOSITION 0
#define DMCP_FIELD_DECOMPOSITION 1
#define DMCP_DOMAIN_DECOMPOSITION 2

#define DMCP_NO_EMBEDDING 0
#define DMCP_FIELD_EMBEDDING 1
#define DMCP_DOMAIN_EMBEDDING 2

struct DM_CP
{
  System *sys;
  std::map<std::string, unsigned int> *varids;
  std::map<unsigned int, std::string> *varnames;
  std::map<std::string, unsigned int> *blockids;
  std::map<unsigned int, std::string> *blocknames;
  unsigned int decomposition_type;
  std::vector<std::set<unsigned int>> *decomposition;
  unsigned int embedding_type;
  IS embedding;
  unsigned int vec_count;
};

struct DMVec_CP
{
  std::string label;
};

#undef __FUNCT__
#define __FUNCT__ "DMCPGetVec_Private"
PetscErrorCode DMCPGetVec_Private(DM /*dm*/, const char * /*name*/, Vec *)
{
  PetscFunctionBegin;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode DMCPSetUpName_Private(DM dm);

// this function is to set dlm->blockids dlm->blocknames;
#undef __FUNCT__
#define __FUNCT__ "DMCPSetSystem_CP"
PetscErrorCode DMCPSetSystem_CP(DM dm, System &sys)
{
  const communicator &comm = sys.comm();

  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  PetscBool iscp;
  ierr = PetscObjectTypeCompare((PetscObject)dm, DMCP, &iscp);
  CHKERRQ(ierr);
  if (!iscp)
    SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Got DM of type %s, not of type %s", ((PetscObject)dm)->type_name, DMCP);

  if (dm->setupcalled)
    SETERRQ(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONGSTATE, "Cannot reset the CP system after DM has been set up.");
  DM_CP *dlm = (DM_CP *)(dm->data);
  dlm->sys = &sys;
  /* Initially populate the sets of active blockids and varids using all of the
     existing blocks/variables (only variables are supported at the moment). */
  DofMap &dofmap = dlm->sys->dof_map();
  dlm->varids->clear();
  dlm->varnames->clear();
  Variables vars = dofmap.variables();
  for (int i = 0; i < vars.size(); i++)
  {
    dlm->varids->insert(std::pair<std::string, unsigned int>(vars.variable_name(i), i));
    dlm->varnames->insert(std::pair<unsigned int, std::string>(i, vars.variable_name(i)));
  }
  const Mesh &mesh = dlm->sys->get_mesh();
  dlm->blockids->clear();
  dlm->blocknames->clear();
  std::set<int> blocks;
  /* The following effectively is a verbatim copy of MeshBase::n_subdomains(). */
  // This requires an inspection on every processor
  // parallel_only(mesh.comm());

  for (const auto elem : mesh.elem_local_range())
  {
    blocks.insert(elem.processor_id());
  }
  // Some subdomains may only live on other processors
  // comm.set_union(blocks);
  std::set<int>::iterator bit = blocks.begin();
  std::set<int>::iterator bend = blocks.end();
  if (bit == bend)
    SETERRQ(((PetscObject)dm)->comm, PETSC_ERR_PLIB, "No mesh blocks found.");
  for (; bit != bend; ++bit)
  {
    int bid = *bit;
    std::string bname = "subdomain_" + std::to_string(bid);
    if (!bname.length())
    {

      std::ostringstream ss;
      ss << "dm" << bid;
      bname = ss.str();
    }
    dlm->blockids->insert(std::pair<std::string, unsigned int>(bname, bid));
    dlm->blocknames->insert(std::pair<unsigned int, std::string>(bid, bname));
  }
  ierr = DMCPSetUpName_Private(dm);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMCPGetSystem_CP"
PetscErrorCode DMCPGetSystem_CP(DM dm, System *&sys)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  PetscBool iscp;
  ierr = PetscObjectTypeCompare((PetscObject)dm, DMCP, &iscp);
  CHKERRQ(ierr);
  if (!iscp)
    SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Got DM of type %s, not of type %s", ((PetscObject)dm)->type_name, DMCP);
  DM_CP *dlm = (DM_CP *)(dm->data);
  sys = dlm->sys;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMCPGetBlocks"
PetscErrorCode DMCPGetBlocks(DM dm, PetscInt *n, char ***blocknames)
{
  PetscErrorCode ierr;
  PetscInt i;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  PetscBool iscp;
  ierr = PetscObjectTypeCompare((PetscObject)dm, DMCP, &iscp);
  CHKERRQ(ierr);
  if (!iscp)
    SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Got DM of type %s, not of type %s", ((PetscObject)dm)->type_name, DMCP);
  DM_CP *dlm = (DM_CP *)(dm->data);
  PetscValidPointer(n, 2);
  *n = static_cast<unsigned int>(dlm->blockids->size());
  if (!blocknames)
    PetscFunctionReturn(static_cast<PetscErrorCode>(0));
  ierr = PetscMalloc(*n * sizeof(char *), blocknames);
  CHKERRQ(ierr);
  i = 0;
  for (const auto &pr : *(dlm->blockids))
  {
    ierr = PetscStrallocpy(pr.first.c_str(), *blocknames + i);
    CHKERRQ(ierr);
    ++i;
  }
  PetscFunctionReturn(static_cast<PetscErrorCode>(0));
}

#undef __FUNCT__
#define __FUNCT__ "DMCPGetVariables"
PetscErrorCode DMCPGetVariables(DM dm, PetscInt *n, char ***varnames)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  PetscBool iscp;
  PetscInt i;
  ierr = PetscObjectTypeCompare((PetscObject)dm, DMCP, &iscp);
  CHKERRQ(ierr);
  if (!iscp)
    SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Got DM of type %s, not of type %s", ((PetscObject)dm)->type_name, DMCP);
  DM_CP *dlm = (DM_CP *)(dm->data);

  *n = static_cast<unsigned int>(dlm->varids->size());
  if (!varnames)
    PetscFunctionReturn(static_cast<PetscErrorCode>(0));
  ierr = PetscMalloc(*n * sizeof(char *), varnames);
  CHKERRQ(ierr);
  i = 0;
  for (const auto &pr : *(dlm->varids))
  {
    ierr = PetscStrallocpy(pr.first.c_str(), *varnames + i);
    CHKERRQ(ierr);
    ++i;
  }
  PetscFunctionReturn(static_cast<PetscErrorCode>(0));
}

// set the block name and the field name;
#undef __FUNCT__
#define __FUNCT__ "DMCPSetUpName_Private"
PetscErrorCode DMCPSetUpName_Private(DM dm)
{
  DM_CP *dlm = (DM_CP *)dm->data;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  std::string name = dlm->sys->name();
  std::map<unsigned int, std::string> *dnames = PETSC_NULLPTR,
                                      *enames = PETSC_NULLPTR;
  if (dlm->decomposition_type == DMCP_FIELD_DECOMPOSITION)
  {
    name += ":dec:var:";
    dnames = dlm->varnames;
  }
  if (dlm->decomposition_type == DMCP_DOMAIN_DECOMPOSITION)
  {
    name += ":dec:block:";
    dnames = dlm->blocknames;
  }
  if (dnames)
  {
    for (auto decomp : *dlm->decomposition)
    {
      for (std::set<unsigned int>::iterator dit_begin = decomp.begin(),
                                            dit = dit_begin,
                                            dit_end = decomp.end();
           dit != dit_end; ++dit)
      {
        unsigned int id = *dit;
        if (dit != dit_begin)
          name += ",";
        name += (*dnames)[id];
      }
      name += ";";
    }
  }
  if (dlm->embedding_type == DMCP_FIELD_EMBEDDING)
  {
    name += ":emb:var:";
    enames = dlm->varnames;
  }
  if (dlm->embedding_type == DMCP_DOMAIN_EMBEDDING)
  {
    name += ":emb:block:";
    enames = dlm->blocknames;
  }
  if (enames)
  {
    for (auto eit = enames->begin(),
              eit_end = enames->end();
         eit != eit_end; ++eit)
    {
      std::string &ename = eit->second;
      if (eit != enames->begin())
        name += ",";
      name += ename;
    }
    name += ";";
  }
  ierr = PetscObjectSetName((PetscObject)dm, name.c_str());
  CHKERRQ(ierr);
  PetscFunctionReturn(static_cast<PetscErrorCode>(0));
}

#undef __FUNCT__
#define __FUNCT__ "DMCreateFieldDecomposition_CP"
static PetscErrorCode DMCreateFieldDecomposition_CP(DM dm, PetscInt *len, char ***namelist, IS **islist, DM **dmlist)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  DM_CP *dlm = (DM_CP *)(dm->data);
  System *sys = dlm->sys;
  Mesh &mesh = sys->get_mesh();
  IS emb;
  if (dlm->decomposition_type != DMCP_FIELD_DECOMPOSITION)
    PetscFunctionReturn(0);

  *len = static_cast<unsigned int>(dlm->decomposition->size());
  if (namelist)
  {
    ierr = PetscMalloc(*len * sizeof(char *), namelist);
    CHKERRQ(ierr);
  }
  if (islist)
  {
    ierr = PetscMalloc(*len * sizeof(IS), islist);
    CHKERRQ(ierr);
  }
  if (dmlist)
  {
    ierr = PetscMalloc(*len * sizeof(DM), dmlist);
    CHKERRQ(ierr);
  }
  DofMap &dofmap = dlm->sys->dof_map();
  for (int i = 0; i < dlm->decomposition->size(); i++)
  {
    std::set<int> dindices;
    std::string dname;
    std::map<std::string, unsigned int> dvarids;
    std::map<unsigned int, std::string> dvarnames;
    unsigned int dvcount = 0;
    for (const auto &v : (*dlm->decomposition)[i])
    {
      std::string vname = (*dlm->varnames)[v];
      dvarids.insert(std::pair<std::string, unsigned int>(vname, v));
      dvarnames.insert(std::pair<unsigned int, std::string>(v, vname));
      if (!dvcount)
        dname = vname;
      else
        dname += "_" + vname;
      ++dvcount;
      if (!islist)
        continue;
      // Iterate only over this DM's blocks.
      for (const auto &pr : *(dlm->blockids))
      {
        const int sbd_id = static_cast<int>(pr.second);
        for (const auto elem : mesh.elem_local_range())
        {
          std::vector<int> evindices;
          // Get the degree of freedom indices for the given variable off the current element.
          int dof = dofmap.dof_indices(elem, v);
          if (dof >= dofmap.first_dof() && dof < dofmap.end_dof()) // might want to use variable_first/last_local_dof instead
            dindices.insert(dof);
        }
      }
    }
    if (namelist)
    {
      ierr = PetscStrallocpy(dname.c_str(), (*namelist) + i);
      CHKERRQ(ierr);
    }
    if (islist)
    {
      IS dis;
      PetscInt *darray;
      ierr = PetscMalloc(sizeof(PetscInt) * dindices.size(), &darray);
      CHKERRQ(ierr);
      int i = 0;
      for (const auto &id : dindices)
      {
        darray[i] = id;
        ++i;
      }
      ierr = ISCreateGeneral(((PetscObject)dm)->comm,
                             static_cast<PetscInt>(dindices.size()),
                             darray, PETSC_OWN_POINTER, &dis);
      CHKERRQ(ierr);
      if (dlm->embedding)
      {
        /* Create a relative embedding into the parent's index space. */
        ierr = ISEmbed(dis, dlm->embedding, PETSC_TRUE, &emb);
        CHKERRQ(ierr);
        PetscInt elen, dlen;
        ierr = ISGetLocalSize(emb, &elen);
        CHKERRQ(ierr);
        ierr = ISGetLocalSize(dis, &dlen);
        CHKERRQ(ierr);
        if (elen != dlen)
          SETERRQ1(((PetscObject)dm)->comm, PETSC_ERR_PLIB, "Failed to embed subdomain %zu", i);
        ierr = ISDestroy(&dis);
        CHKERRQ(ierr);
        dis = emb;
      }
      else
      {
        emb = dis;
      }
      (*islist)[i] = dis;
    }
    if (dmlist)
    {
      DM ddm;
      ierr = DMCreate(((PetscObject)dm)->comm, &ddm);
      CHKERRQ(ierr);
      ierr = DMSetType(ddm, DMCP);
      CHKERRQ(ierr);
      DM_CP *ddlm = (DM_CP *)(ddm->data);
      ddlm->sys = dlm->sys;
      /* copy over the block ids and names */
      *ddlm->blockids = *dlm->blockids;
      *ddlm->blocknames = *dlm->blocknames;
      /* set the vars from the d-th part of the decomposition. */
      *ddlm->varids = dvarids;
      *ddlm->varnames = dvarnames;
      ierr = PetscObjectReference((PetscObject)emb);
      CHKERRQ(ierr);
      ddlm->embedding = emb;
      ddlm->embedding_type = DMCP_FIELD_EMBEDDING;

      ierr = DMCPSetUpName_Private(ddm);
      CHKERRQ(ierr);
      ierr = DMSetFromOptions(ddm);
      CHKERRQ(ierr);
      (*dmlist)[i] = ddm;
    }
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMCreateDomainDecomposition_CP"
static PetscErrorCode DMCreateDomainDecomposition_CP(DM dm, PetscInt *len, char ***namelist, IS **innerislist, IS **outerislist, DM **dmlist)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  DM_CP *dlm = (DM_CP *)(dm->data);
  System *sys = dlm->sys;
  IS emb;
  if (dlm->decomposition_type != DMCP_DOMAIN_DECOMPOSITION)
    PetscFunctionReturn(0);
  *len = static_cast<unsigned int>(dlm->decomposition->size());
  if (namelist)
  {
    ierr = PetscMalloc(*len * sizeof(char *), namelist);
    CHKERRQ(ierr);
  }
  if (innerislist)
  {
    ierr = PetscMalloc(*len * sizeof(IS), innerislist);
    CHKERRQ(ierr);
  }
  if (outerislist)
    *outerislist = PETSC_NULLPTR; /* FIX: allow mesh-based overlap. */
  if (dmlist)
  {
    ierr = PetscMalloc(*len * sizeof(DM), dmlist);
    CHKERRQ(ierr);
  }
  for (int i = 0; i < (*dlm).decomposition->size(); i++)
  {
    std::set<unsigned int> dindices;
    std::string dname;
    std::map<std::string, unsigned int> dblockids;
    std::map<unsigned int, std::string> dblocknames;
    unsigned int dbcount = 0;
    for (const auto &b : (*dlm->decomposition)[i])
    {
      std::string bname = (*dlm->blocknames)[b];
      dblockids.insert(std::pair<std::string, unsigned int>(bname, b));
      dblocknames.insert(std::pair<unsigned int, std::string>(b, bname));
      if (!dbcount)
        dname = bname;
      else
        dname += "_" + bname;
      ++dbcount;
      if (!innerislist)
        continue;
      const PetscInt b_sbd_id = static_cast<PetscInt>(b);
      // Iterate only over this DM's variables.
      // Get the degree of freedom indices for the given variable off the current element.
      for (const auto &elem : sys->get_mesh().elem_local_range())
      {
        int evindices;
        for (const auto &pr : *(dlm->varids))
        {
          // Get the degree of freedom indices for the given variable off the current element.
          sys->dof_map().dof_indices(elem, evindices);

          if (evindices >= sys->dof_map().first_dof() && evindices < sys->dof_map().end_dof()) // might want to use variable_first/last_local_dof instead
            dindices.insert(evindices);
        }
      }
    }
    if (namelist)
    {
      ierr = PetscStrallocpy(dname.c_str(), (*namelist) + i);
      CHKERRQ(ierr);
    }
    if (innerislist)
    {
      PetscInt *darray;
      IS dis;
      ierr = PetscMalloc(sizeof(PetscInt) * dindices.size(), &darray);
      CHKERRQ(ierr);
      int i = 0;
      for (const auto &id : dindices)
      {
        darray[i] = id;
        ++i;
      }
      ierr = ISCreateGeneral(((PetscObject)dm)->comm,
                             static_cast<PetscInt>(dindices.size()),
                             darray, PETSC_OWN_POINTER, &dis);
      CHKERRQ(ierr);
      if (dlm->embedding)
      {
        /* Create a relative embedding into the parent's index space. */
        ierr = ISEmbed(dis, dlm->embedding, PETSC_TRUE, &emb);
        CHKERRQ(ierr);
        PetscInt elen, dlen;
        ierr = ISGetLocalSize(emb, &elen);
        CHKERRQ(ierr);
        ierr = ISGetLocalSize(dis, &dlen);
        CHKERRQ(ierr);
        if (elen != dlen)
          SETERRQ1(((PetscObject)dm)->comm, PETSC_ERR_PLIB, "Failed to embed field %zu", i);
        ierr = ISDestroy(&dis);
        CHKERRQ(ierr);
        dis = emb;
      }
      else
      {
        emb = dis;
      }
      if (innerislist)
      {
        ierr = PetscObjectReference((PetscObject)dis);
        CHKERRQ(ierr);
        (*innerislist)[i] = dis;
      }
      ierr = ISDestroy(&dis);
      CHKERRQ(ierr);
    }
    if (dmlist)
    {
      DM ddm;
      ierr = DMCreate(((PetscObject)dm)->comm, &ddm);
      CHKERRQ(ierr);
      ierr = DMSetType(ddm, DMCP);
      CHKERRQ(ierr);
      DM_CP *ddlm = (DM_CP *)(ddm->data);
      ddlm->sys = dlm->sys;
      /* copy over the varids and varnames */
      *ddlm->varids = *dlm->varids;
      *ddlm->varnames = *dlm->varnames;
      /* set the blocks from the d-th part of the decomposition. */
      *ddlm->blockids = dblockids;
      *ddlm->blocknames = dblocknames;
      ierr = PetscObjectReference((PetscObject)emb);
      CHKERRQ(ierr);
      ddlm->embedding = emb;
      ddlm->embedding_type = DMCP_DOMAIN_EMBEDDING;

      ierr = DMCPSetUpName_Private(ddm);
      CHKERRQ(ierr);
      ierr = DMSetFromOptions(ddm);
      CHKERRQ(ierr);
      (*dmlist)[i] = ddm;
    }
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMCPCreateFieldDecompositionDM"
PetscErrorCode DMCPCreateFieldDecompositionDM(DM dm, PetscInt dnumber, PetscInt *dsizes, char ***dvarlists, DM *ddm)
{
  PetscErrorCode ierr;
  PetscBool isCP;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  ierr = PetscObjectTypeCompare((PetscObject)dm, DMCP, &isCP);
  CHKERRQ(ierr);
  if (!isCP)
    SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Got DM of type %s, not of type %s", ((PetscObject)dm)->type_name, DMCP);
  if (dnumber < 0)
    SETERRQ1(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Negative number %" PetscInt_FMT " of decomposition parts", dnumber);
  // #if PETSC_RELEASE_GREATER_EQUALS(3,20,0)
  //   PetscAssertPointer(ddm,5);

  PetscValidPointer(ddm, 5);
  DM_CP *dlm = (DM_CP *)(dm->data);
  ierr = DMCreate(((PetscObject)dm)->comm, ddm);
  CHKERRQ(ierr);
  ierr = DMSetType(*ddm, DMCP);
  CHKERRQ(ierr);
  DM_CP *ddlm = (DM_CP *)((*ddm)->data);
  ddlm->sys = dlm->sys;
  ddlm->varids = dlm->varids;
  ddlm->varnames = dlm->varnames;
  ddlm->blockids = dlm->blockids;
  ddlm->blocknames = dlm->blocknames;
  ddlm->decomposition = new (std::vector<std::set<unsigned int>>);
  ddlm->decomposition_type = DMCP_FIELD_DECOMPOSITION;
  if (dnumber)
  {
    for (PetscInt d = 0; d < dnumber; ++d)
    {
      if (dsizes[d] < 0)
        SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Negative size %" PetscInt_FMT " of decomposition part %" PetscInt_FMT, dsizes[d], d);
      ddlm->decomposition->push_back(std::set<unsigned int>());
      for (PetscInt v = 0; v < dsizes[d]; ++v)
      {
        std::string vname(dvarlists[d][v]);
        std::map<std::string, unsigned int>::const_iterator vit = dlm->varids->find(vname);
        if (vit == dlm->varids->end())
          SETERRQ3(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Variable %" PetscInt_FMT " on the %" PetscInt_FMT "-th list with name %s is not owned by this DM", v, d, dvarlists[d][v]);
        unsigned int vid = vit->second;
        (*ddlm->decomposition)[d].insert(vid);
      }
    }
  }
  else
  { // Empty splits indicate default: split all variables with one per split.
    PetscInt d = 0;
    for (const auto &pr : (*ddlm->varids))
    {
      ddlm->decomposition->push_back(std::set<unsigned int>());
      unsigned int vid = pr.second;
      (*ddlm->decomposition)[d].insert(vid);
      ++d;
    }
  }
  ierr = DMCPSetUpName_Private(*ddm);
  CHKERRQ(ierr);
  ierr = DMSetFromOptions(*ddm);
  CHKERRQ(ierr);
  ierr = DMSetUp(*ddm);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMCPCreateDomainDecompositionDM"
PetscErrorCode DMCPCreateDomainDecompositionDM(DM dm, PetscInt dnumber, PetscInt *dsizes, char ***dblocklists, DM *ddm)
{
  PetscErrorCode ierr;
  PetscBool isCP;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  ierr = PetscObjectTypeCompare((PetscObject)dm, DMCP, &isCP);
  CHKERRQ(ierr);
  if (!isCP)
    SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Got DM of type %s, not of type %s", ((PetscObject)dm)->type_name, DMCP);
  if (dnumber < 0)
    SETERRQ1(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Negative number %"
                                                           " of decomposition parts",
             dnumber);
  // #if PETSC_RELEASE_GREATER_EQUALS(3,20,0)
  //  PetscAssertPointer(ddm,5);

  PetscValidPointer(ddm, 5);

  DM_CP *dlm = (DM_CP *)(dm->data);
  ierr = DMCreate(((PetscObject)dm)->comm, ddm);
  CHKERRQ(ierr);
  ierr = DMSetType(*ddm, DMCP);
  CHKERRQ(ierr);
  DM_CP *ddlm = (DM_CP *)((*ddm)->data);
  ddlm->sys = dlm->sys;
  ddlm->varids = dlm->varids;
  ddlm->varnames = dlm->varnames;
  ddlm->blockids = dlm->blockids;
  ddlm->blocknames = dlm->blocknames;
  ddlm->decomposition = new (std::vector<std::set<unsigned int>>);
  ddlm->decomposition_type = DMCP_DOMAIN_DECOMPOSITION;
  if (dnumber)
  {
    for (PetscInt d = 0; d < dnumber; ++d)
    {
      if (dsizes[d] < 0)
        SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Negative size %" PetscInt_FMT " of decomposition part %" PetscInt_FMT, dsizes[d], d);
      ddlm->decomposition->push_back(std::set<unsigned int>());
      for (PetscInt b = 0; b < dsizes[d]; ++b)
      {
        std::string bname(dblocklists[d][b]);
        std::map<std::string, unsigned int>::const_iterator bit = dlm->blockids->find(bname);
        if (bit == dlm->blockids->end())
          SETERRQ3(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Block %" PetscInt_FMT " on the %" PetscInt_FMT "-th list with name %s is not owned by this DM", b, d, dblocklists[d][b]);
        unsigned int bid = bit->second;
        (*ddlm->decomposition)[d].insert(bid);
      }
    }
  }
  else
  { // Empty splits indicate default: split all blocks with one per split.
    PetscInt d = 0;
    for (const auto &pr : (*ddlm->blockids))
    {
      ddlm->decomposition->push_back(std::set<unsigned int>());
      unsigned int bid = pr.second;
      (*ddlm->decomposition)[d].insert(bid);
      ++d;
    }
  }
  ierr = DMCPSetUpName_Private(*ddm);
  CHKERRQ(ierr);
  ierr = DMSetFromOptions(*ddm);
  CHKERRQ(ierr);
  ierr = DMSetUp(*ddm);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMCreateGlobalVector_CP"
static PetscErrorCode DMCreateGlobalVector_CP(DM dm, Vec *x)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  DM_CP *dlm = (DM_CP *)(dm->data);
  PetscBool eq;

  ierr = PetscObjectTypeCompare((PetscObject)dm, DMCP, &eq);
  CHKERRQ(ierr);

  if (!eq)
    SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "DM of type %s, not of type %s", ((PetscObject)dm)->type_name, DMCP);

  if (!dlm->sys)
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "No cp system set for DM_CP");

  if (dlm->embedding)
  {
    PetscInt n;
    ierr = VecCreate(dlm->sys->comm(), x);
    CHKERRQ(ierr);
    ierr = ISGetLocalSize(dlm->embedding, &n);
    CHKERRQ(ierr);
    ierr = VecSetSizes(*x, n, PETSC_DETERMINE);
    CHKERRQ(ierr);
    ierr = VecSetType(*x, VECMPI);
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(*x);
    CHKERRQ(ierr);
    ierr = VecSetUp(*x);
    CHKERRQ(ierr);
  }
  else
  {
    ierr = dlm->sys->create_global_vec(x);
    CHKERRQ(ierr);
  }
  ierr = VecSetDM(*x, dm);
  CHKERRQ(ierr);
  PetscFunctionReturn(PETSC_SUCCESS);
}

#undef __FUNCT__
#define __FUNCT__ "DMCreateLocalVector_CP"
static PetscErrorCode DMCreateLocalVector_CP(DM dm, Vec *x)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  DM_CP *dlm = (DM_CP *)(dm->data);
  PetscBool eq;

  ierr = PetscObjectTypeCompare((PetscObject)dm, DMCP, &eq);
  CHKERRQ(ierr);

  if (!eq)
    SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "DM of type %s, not of type %s", ((PetscObject)dm)->type_name, DMCP);

  if (!dlm->sys)
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "No cp system set for DM_CP");
  if (dlm->embedding)
  {
    PetscInt n;
    ierr = VecCreate(dlm->sys->comm(), x);
    CHKERRQ(ierr);
    ierr = ISGetLocalSize(dlm->embedding, &n);
    CHKERRQ(ierr);
    ierr = VecSetSizes(*x, n, PETSC_DETERMINE);
    CHKERRQ(ierr);
    ierr = VecSetType(*x, VECMPI);
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(*x);
    CHKERRQ(ierr);
    ierr = VecSetUp(*x);
    CHKERRQ(ierr);
  }
  else
  {
    ierr = dlm->sys->create_local_vec(x);
    CHKERRQ(ierr);
  }
  ierr = VecSetDM(*x, dm);
  CHKERRQ(ierr);
  PetscFunctionReturn(PETSC_SUCCESS);
}

#undef __FUNCT__
#define __FUNCT__ "DMCreateMatrix_CP"
static PetscErrorCode DMCreateMatrix_CP(DM dm, Mat *A)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  DM_CP *dlm = (DM_CP *)(dm->data);
  PetscBool eq;
  ierr = PetscObjectTypeCompare((PetscObject)dm, DMCP, &eq);
  CHKERRQ(ierr);

  if (!eq)
    SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "DM of type %s, not of type %s", ((PetscObject)dm)->type_name, DMCP);

  dlm->sys->create_mat(*A);
  ierr = PetscObjectReference((PetscObject)(*A));
  CHKERRQ(ierr);
  PetscFunctionReturn(static_cast<PetscErrorCode>(0));
}

#undef __FUNCT__
#define __FUNCT__ "DMView_CP"
static PetscErrorCode DMView_CP(DM dm, PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscBool isascii;
  const char *name, *prefix;
  DM_CP *dlm = (DM_CP *)dm->data;
  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer, PETSCVIEWERASCII, &isascii);
  CHKERRQ(ierr);
  if (isascii)
  {
    ierr = PetscObjectGetName((PetscObject)dm, &name);
    CHKERRQ(ierr);
    ierr = PetscObjectGetOptionsPrefix((PetscObject)dm, &prefix);
    CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer, "DM cp with name %s and prefix %s\n", name, prefix);
    CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer, "blocks:");
    CHKERRQ(ierr);
    std::map<std::string, unsigned int>::iterator bit = dlm->blockids->begin();
    std::map<std::string, unsigned int>::const_iterator bend = dlm->blockids->end();
    for (; bit != bend; ++bit)
    {
      ierr = PetscViewerASCIIPrintf(viewer, "(%s,%d) ", bit->first.c_str(), bit->second);
      CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer, "\n");
    CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer, "variables:");
    CHKERRQ(ierr);
    std::map<std::string, unsigned int>::iterator vit = dlm->varids->begin();
    std::map<std::string, unsigned int>::const_iterator vend = dlm->varids->end();
    for (; vit != vend; ++vit)
    {
      ierr = PetscViewerASCIIPrintf(viewer, "(%s,%d) ", vit->first.c_str(), vit->second);
      CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer, "\n");
    CHKERRQ(ierr);
    if (dlm->decomposition_type == 0)
    {
      ierr = PetscViewerASCIIPrintf(viewer, "No decomposition\n");
      CHKERRQ(ierr);
    }
    else
    {
      if (dlm->decomposition_type == DMCP_FIELD_DECOMPOSITION)
      {
        ierr = PetscViewerASCIIPrintf(viewer, "Field decomposition by variable: ");
        CHKERRQ(ierr);
      }
      else if (dlm->decomposition_type == DMCP_DOMAIN_DECOMPOSITION)
      {
        ierr = PetscViewerASCIIPrintf(viewer, "Domain decomposition by block: ");
        CHKERRQ(ierr);
      }
      else
        SETERRQ1(((PetscObject)dm)->comm, PETSC_ERR_PLIB, "Unexpected decomposition type: %d", dlm->decomposition_type);
      /* FIX: decompositions might have different sizes and components on different ranks. */
      for (auto &set : *dlm->decomposition)
      {
        auto it = set.begin();
        auto end = set.end();
        if (it != end)
        {
          ierr = PetscViewerASCIIPrintf(viewer, "%u", *it);
          CHKERRQ(ierr);
          ++it;
        }
        for (; it != end; ++it)
        {
          ierr = PetscViewerASCIIPrintf(viewer, ",%u", *it);
          CHKERRQ(ierr);
        }
        ierr = PetscViewerASCIIPrintf(viewer, ";");
        CHKERRQ(ierr);
      }

      ierr = PetscViewerASCIIPrintf(viewer, "\n");
      CHKERRQ(ierr);
    }
  }

  PetscFunctionReturn(static_cast<PetscErrorCode>(0));
}

#undef __FUNCT__
#define __FUNCT__ "DMSetUp_libMesh"
static PetscErrorCode DMSetUp_CP(DM dm)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  DM_CP *dlm = (DM_CP *)(dm->data);
  PetscBool eq;

  ierr = PetscObjectTypeCompare((PetscObject)dm, DMCP, &eq);
  CHKERRQ(ierr);

  if (!eq)
    SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "DM of type %s, not of type %s", ((PetscObject)dm)->type_name, DMCP);

  if (!dlm->sys)
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "No corner_point system set for DM_CP");

  if (dlm->embedding)
  {
    dm->ops->createglobalvector = 0;
    dm->ops->creatematrix = 0;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDestroy_libMesh"
static PetscErrorCode DMDestroy_CP(DM dm)
{
  DM_CP *dlm = (DM_CP *)(dm->data);
  PetscErrorCode ierr;
  PetscFunctionBegin;
  delete dlm->varids;
  delete dlm->varnames;
  delete dlm->blockids;
  delete dlm->blocknames;
  delete dlm->decomposition;
  ierr = ISDestroy(&dlm->embedding);
  CHKERRQ(ierr);
  ierr = PetscFree(dm->data);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode DMCreate_CP(DM dm)
{
  PetscErrorCode ierr;
  DM_CP *dlm;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  ierr = PetscNew(&dlm);
  CHKERRQ(ierr);
  dm->data = dlm;

  /* DMCP impl */
  dlm->varids = new (std::map<std::string, unsigned int>);
  dlm->blockids = new (std::map<std::string, unsigned int>);
  dlm->varnames = new (std::map<unsigned int, std::string>);
  dlm->blocknames = new (std::map<unsigned int, std::string>);
  dlm->decomposition = PETSC_NULLPTR;
  dlm->decomposition_type = DMCP_NO_DECOMPOSITION;

  /* DM API */

  dm->ops->createglobalvector = DMCreateGlobalVector_CP;
  dm->ops->createlocalvector = DMCreateLocalVector_CP;
  dm->ops->getcoloring = 0;
  dm->ops->creatematrix = DMCreateMatrix_CP;
  dm->ops->createinterpolation = 0; //
  dm->ops->refine = 0;
  dm->ops->coarsen = 0;
  // * dm->ops->getinjection was renamed to dm->ops->createinjection in PETSc 5a84ad338 (5 Jul 2019)
  // * dm->ops-getaggregates was removed in PETSc 97779f9a (5 Jul 2019)
  // * Both changes were merged into PETSc master in 94aad3ce (7 Jul 2019).
  PetscInt major, minor, subminor, release;
  PetscCall(PetscGetVersionNumber(&major, &minor, &subminor, &release));

  dm->ops->createinjection = 0;

  dm->ops->createfielddecomposition = DMCreateFieldDecomposition_CP;
  dm->ops->createdomaindecomposition = DMCreateDomainDecomposition_CP;
  dm->ops->destroy = DMDestroy_CP;
  dm->ops->view = DMView_CP;
  dm->ops->setfromoptions = 0; // DMSetFromOptions_CP;
  dm->ops->setup = DMSetUp_CP;
  /* DMCP API */
  ierr = PetscObjectComposeFunction((PetscObject)dm, "DMCPSetSystem_C", DMCPSetSystem_CP);
  CHKERRQ(ierr);

  ierr = PetscObjectComposeFunction((PetscObject)dm, "DMCPGetSystem_C", DMCPGetSystem_CP);
  CHKERRQ(ierr);
  PetscFunctionReturn(static_cast<PetscErrorCode>(0));
}

EXTERN_C_END
