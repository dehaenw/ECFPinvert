import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

#atoms dict and valence
ATOMSDICT = {6:4,7:3,8:2,9:1} #CNOF
ATOMSDICT_CATION = {7:4,8:3}
ATOMSDICT_ANION = {7:2,8:1}
ATOMSDICT_FULL = {6:4,7:3,8:2,9:1,15:5,16:6,17:1,35:1,53:1} #CNOFPSClBrI
BONDTYPES = {"s":Chem.BondType.SINGLE,"d":Chem.BondType.DOUBLE,"t":Chem.BondType.TRIPLE,"a":Chem.BondType.AROMATIC}


def gen_radius_0_molfrags(atomsdict):
    """ generate all fragments based on the atoms defined in the atoms dict
    """
    frags = []
    for anumber in atomsdict.keys():
        #acyclic frags:
        pb = get_possible_bonds(atomsdict[anumber])
        for curr_bonds in pb: 
            mol = Chem.RWMol(Chem.Mol())
            mol.AddAtom(Chem.Atom(anumber)) #this returns the RDKit index of the atom
            for curr_bond in curr_bonds:
                idx2 = mol.AddAtom(Chem.Atom(0))
                mol.GetAtomWithIdx(idx2).SetProp('cyclic','no')
                bondIdx = mol.AddBond(0,idx2, BONDTYPES[curr_bond])
            frags.append(mol)
        #cyclic frags
        pb = get_possible_bonds(atomsdict[anumber],inring=True)
        for curr_bonds in pb: 
            mol = Chem.RWMol(Chem.Mol())
            mol.AddAtom(Chem.Atom(anumber)) #this returns the RDKit index of the atom
            for i,curr_bond in enumerate(curr_bonds):
                idx2 = mol.AddAtom(Chem.Atom(0))
                if i<2 or curr_bond=="a": #aromatic must always be in ring
                    mol.GetAtomWithIdx(idx2).SetProp('cyclic','yes') #first two dummies must be in ring
                else:
                    mol.GetAtomWithIdx(idx2).SetProp('cyclic','maybe') #next ones can be, but dont have to be in ring
                mol.AddBond(0,idx2, BONDTYPES[curr_bond])
            mol.AddBond(1,2, Chem.BondType.SINGLE)
            if len(curr_bonds) >= 3:
                mol.AddBond(2,3, Chem.BondType.SINGLE)
            frags.append(mol)
            if anumber in [7,15,33,51] and curr_bonds == ["a","a"]: # need to explicitly put [nH],[pH] etc bc n1cccc1 eg doesnt parts needs to be [nH]1cccc1
                extrafrag = Chem.RWMol(mol)
                extrafrag.GetAtomWithIdx(0).SetNumExplicitHs(1)
                frags.append(extrafrag)
    #this operation has to be done because during sanitization some reassignment of bond type (eg between aromatic and single) occutrs due to how duimmies are handled
    frags = [sanitize_aro(fr) for fr in frags]
    return frags
    
def gen_charged_frags():
    """ add charged fragments. this is because some charges species are quite
    common in sets such as chembl. the charged species of interest are: N+ N-
    O- and O+ with various bondtypes hcounts etc
    """
    frags = []
    for j,atomsdict in enumerate([ATOMSDICT_ANION,ATOMSDICT_CATION]):
        for anumber in atomsdict.keys():
            #acyclic frags:
            pb = get_possible_bonds(atomsdict[anumber])
            for curr_bonds in pb: 
                mol = Chem.RWMol(Chem.Mol())
                mol.AddAtom(Chem.Atom(anumber)) #this returns the RDKit index of the atom
                mol.GetAtoms()[0]
                mol.GetAtoms()[0].SetFormalCharge([-1,1][j])
                for curr_bond in curr_bonds:
                    idx2 = mol.AddAtom(Chem.Atom(0))
                    mol.GetAtomWithIdx(idx2).SetProp('cyclic','no')
                    bondIdx = mol.AddBond(0,idx2, BONDTYPES[curr_bond])
                frags.append(mol)
            #cyclic frags
            pb = get_possible_bonds(atomsdict[anumber],inring=True)
            for curr_bonds in pb: 
                mol = Chem.RWMol(Chem.Mol())
                mol.AddAtom(Chem.Atom(anumber)) #this returns the RDKit index of the atom
                mol.GetAtoms()[0]
                mol.GetAtoms()[0].SetFormalCharge([-1,1][j])
                for i,curr_bond in enumerate(curr_bonds):
                    idx2 = mol.AddAtom(Chem.Atom(0))
                    if i<2 or curr_bond=="a": #aromatic must always be in ring
                        mol.GetAtomWithIdx(idx2).SetProp('cyclic','yes') #first two dummies must be in ring
                    else:
                        mol.GetAtomWithIdx(idx2).SetProp('cyclic','maybe') #next ones can be, but dont have to be in ring
                    mol.AddBond(0,idx2, BONDTYPES[curr_bond])
                mol.AddBond(1,2, Chem.BondType.SINGLE)
                if len(curr_bonds) >= 3:
                    mol.AddBond(2,3, Chem.BondType.SINGLE)
                frags.append(mol)
        #this operation has to be done because during sanitization some reassignment of bond type (eg between aromatic and single) occutrs due to how duimmies are handled
    charged_frags = [sanitize_aro(fr) for fr in frags if sanitize_aro(fr)!=None]
    return charged_frags

def sanitize_aro(mol):
    """ rdkit automatically converts stuff like *1O*1 to aromatic when
    sanitizingm this function is a way to subvert this.
    """
    try:
        Chem.SanitizeMol(mol,sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_KEKULIZE^Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
    except:
        print("ERROR !!! sanitization without aromatic changes did not work!!!! this is the aromaticity for each atom")
        for at in mol.GetAtoms():
            print(at.GetIsAromatic())
        print("and the bonds")
        for bo in mol.GetBonds():
            print(bo.GetBondType())
        print("and the problem smiles")
        print(Chem.MolToSmiles(mol))
    return mol


def inspect_frag(frag):
    """get some useful info out of a radius 0 fragment:
    """
    fraginfo = {}
    fraginfo["ringcount"]=frag.GetRingInfo().NumRings()
    fraginfo["valence"] = 0
    fraginfo["bondtypes"] = []
    for at in frag.GetAtoms():
        if at.GetAtomicNum()==0:
            fraginfo["valence"]+=1
            fraginfo["bondtypes"].append(frag.GetBondBetweenAtoms(0,at.GetIdx()).GetBondType())
    return fraginfo

def frag_ringopen(frag,valence,setthird=False):
    """ helper function for when a new ring is opened by adding a valence 3+
    fragment
    """
    frag_prep = Chem.RWMol(frag)
    for a in range(valence):
        for b in range(valence):
            if b>a:
                frag_prep.RemoveBond(a+1, b+1)
    if setthird == True:
        frag_prep.GetAtomWithIdx(3).SetProp('cyclic','yes')
    return frag_prep
        
def get_possible_bonds(valence,inring=False):
    """ helper function to enumerate the possible bond hcount etc options
    for each atom based on its valence
    """
    possible_bonds = []
    if inring==False:
        possible_bonds += [["s"]]
        if valence >= 2:
            possible_bonds += [["s","s"],["d"]]
        if valence >= 3:
            possible_bonds += [["s","s","s"],["d","s"],["s","d"],["t"]]
        if valence >= 4:
            possible_bonds += [["s","s","s","s"],["d","s","s"],["s","d","s"],["d","d"],["t","s"],["s","t"]]
        if valence >= 5:
            possible_bonds += [["s","d","d"],["d","d","s"],["s","s","s","s","s"],["d","s","s","s"],["s","d","s","s"]] #ie phosphorous
        if valence >= 6:
            possible_bonds += [["s","d","d","s"],["d","d","s","s"],["s","s","s","s","s","s"],["d","d","d"]] #ie sulfur
    else:
        if valence >= 2:
            possible_bonds += [["s","s"],["a","a"]] # [a,a] bc of oxygen eg in o1cccc1
        if valence >= 3:
            possible_bonds += [["s","s","s"],["d","s"],["s","d"],["a","a","s"],["s","a","a"],["a","a","a"]] #add ["s","a","a"],["a","s","a"] ?
        if valence >= 4:
            possible_bonds += [["s","s","s","s"],["d","s","s"],["s","d","s"],["s","s","d"],["d","d"],["t","s"],["s","t"],["a","a","d"],["d","a","a"]] #aad and daa added because ofOn1cnc(=N)[nH]1 type systems
        if valence >= 5:
            possible_bonds += [["s","d","d"],["d","d","s"],["d","s","d"],["s","s","s","s","s"],["s","s","d","s"],["s","d","s","s"]]          
        if valence >= 6:
            possible_bonds += [["s","s","d","d"],["d","d","s","s"],["d","s","d","s"],["s","d","d","s"],["d","d","d"]]
    return possible_bonds
    
    
def add_fs_mol_ro(mol,idx,a,alen,frag,infolist):
    """helper function for build_mol_from_fs used to open a new ring.
    """
    if infolist[idx]["valence"]>=3:
        ro_frag = frag_ringopen(frag,infolist[idx]["valence"])
        mol = Chem.RWMol(Chem.CombineMols(mol,ro_frag))
        abondtype = mol.GetBondBetweenAtoms(a-1,a).GetBondType()
        bbondtype = infolist[idx]["bondtypes"][0]
        if bbondtype == abondtype:
            mol.AddBond(a-1,alen,abondtype)
            mol.RemoveAtom(alen+1)
            mol.RemoveAtom(a)
            for k in range(infolist[idx]["valence"]-1):
                mol.AddBond(0,alen+k,infolist[idx]["bondtypes"][k+1])
        else:
            mol=None #bonds must match
    else:
        # if there is more than 1 dummies, open ring between dummy 1 and dummy 2
        dummys = mol.GetSubstructMatches(Chem.MolFromSmarts("[#0]"))
        if  len(dummys)>= 2:
            ro_frag = frag_ringopen(frag,infolist[idx]["valence"])
            mol = Chem.RWMol(Chem.CombineMols(mol,ro_frag))
            abondtype = mol.GetBondBetweenAtoms(a-1,a).GetBondType()
            bbondtype = infolist[idx]["bondtypes"][0]
            if bbondtype == abondtype:
                mol.AddBond(a-1,alen,abondtype)
                mol.AddBond(0,alen+2,infolist[idx]["bondtypes"][1])
                mol.RemoveAtom(alen+1)
                mol.RemoveAtom(a)
            else:
                mol=None #bonds must match
        else:    
            mol=None #need at least 3 bonds, 1 to connect to dummy and 2 for cycle
    return mol
    
def add_fs_mol_rcont(mol,idx,a,alen,frag,infolist):
    """ helper function for build_mol_from_fs used to continue building
    inside an unclosed ring.
    """
    ro_frag = frag_ringopen(frag,infolist[idx]["valence"])
    mol = Chem.RWMol(Chem.CombineMols(mol,ro_frag))
    abondtype = mol.GetBondBetweenAtoms(a-1,a).GetBondType()
    bbondtype = infolist[idx]["bondtypes"][0]
    if bbondtype == abondtype:
        mol.AddBond(a-1,alen,abondtype)
        mol.AddBond(0,alen+2,infolist[idx]["bondtypes"][1])
        for k in range(max(0,infolist[idx]["valence"]-2)):
            mol.AddBond(0,alen+k+3,infolist[idx]["bondtypes"][k+2]) 
            mol.GetAtomWithIdx(alen+k+3).SetProp("cyclic","maybe")
        mol.RemoveAtom(alen+1)
        mol.RemoveAtom(a)
    else:
        mol=None #bonds must match  
    return mol
    
def add_fs_mol_nr(mol,idx,a,alen,frag,infolist):
    """ helper function for build_mol_from_fs used to add a non ring/acyclic
    atom
    """
    if mol.GetAtomWithIdx(a).GetProp("cyclic") != "yes":
        mol = Chem.RWMol(Chem.CombineMols(mol,frag))
        abondtype = mol.GetBondBetweenAtoms(a-1,a).GetBondType()
        bbondtype = infolist[idx]["bondtypes"][0]
        if bbondtype == abondtype:
            mol.AddBond(a-1,alen,abondtype)
            mol.RemoveAtom(alen+1)
            mol.RemoveAtom(a)
        else:
            mol = None #bonds must match 
    else:
        mol = None #dummy is in cycle so cant use acyclic frag
    return mol

def add_fs_mol(mol,idx,fraglist,infolist):
    """ helper function for build_mol_from_fs used to add an atom
    (i.e. non ring closure step) can be expanding ring, acyclic frag addition
    or opening a new ring
    """
    frag = fraglist[idx]
    ssm=mol.GetSubstructMatches(Chem.MolFromSmarts("[#0]"))
    if len(ssm)>0: #have to do this check because it crashes if theres a sequences like C-* C-* C-* where there are no more dummies at the 3rd step
        a = ssm[0][0] 
    else:
        a = 0
    alen = len(mol.GetAtoms())
    if a != 0:
        if infolist[idx]["ringcount"]==0:
            mol=add_fs_mol_nr(mol,idx,a,alen,frag,infolist)
        else:
            if mol.GetAtomWithIdx(a).GetProp("cyclic") in ["no","maybe"]:
                mol = add_fs_mol_ro(mol,idx,a,alen,frag,infolist)
            elif mol.GetAtomWithIdx(a).GetProp("cyclic") == "yes":
                mol = add_fs_mol_rcont(mol,idx,a,alen,frag,infolist)                           
    else:
        mol = None
    return mol
    
def add_fs_mol_rc(mol,idx,fraglist,infolist):
    """helper function for build_mol_from_fs used to perform a ring closure.
    this happens when idx is negative. idx -1 means a ring closure from dummy
    1 to dummy 2. idx -2 is closure from dummy 1 to dummy 3. et cetera
    """
    bonds = mol.GetAtomWithIdx(0).GetBonds()
    temp_mol = Chem.RWMol(mol) #need to copy it because ringsize needs to be checked
    temp_mol.RemoveAtom(0)
    dm = Chem.GetDistanceMatrix(temp_mol)
    if -idx<=(len(bonds)-1):
        if dm[bonds[0].GetEndAtomIdx()-1][bonds[-idx].GetEndAtomIdx()-1]>3: 
            d1idx=bonds[0].GetEndAtomIdx()
            d2idx=bonds[-idx].GetEndAtomIdx() 
            bt1=mol.GetAtomWithIdx(d1idx).GetBonds()[0].GetBondType()
            bt2=mol.GetAtomWithIdx(d2idx).GetBonds()[0].GetBondType()
            nondummyidx=mol.GetAtomWithIdx(d1idx).GetBonds()[0].GetBeginAtomIdx()
            nondummyidx2=mol.GetAtomWithIdx(d2idx).GetBonds()[0].GetBeginAtomIdx()
            if bt1==bt2:
                mol.AddBond(nondummyidx,nondummyidx2,bt1) 
                mol.RemoveBond(nondummyidx2,d2idx)
                mol.RemoveBond(nondummyidx,d1idx)
                mol.RemoveBond(0,d2idx)
                mol.RemoveBond(0,d1idx)
                for idx in sorted([d1idx,d2idx])[::-1]:
                    mol.RemoveAtom(idx)
                for at in mol.GetAtoms():
                    if at.GetAtomicNum()==0:
                        if at.GetProp("cyclic")=="yes":
                            pass# check if it should be replaced with at.SetProp("cyclic","maybe") #otherwise it can get stuck in a bad half ring closure liek CC1(C)CC1(C)O
            else:
                mol=None
        else:
            mol=None
    else:
        mol=None
    return mol

def build_mol_from_fs(fs,fraglist,infolist):
    """ given a list of 1 atom fragment and ring closure operations construct
    a molecule. return None in case of any issues along the way.
    """
    fi=inspect_frag(fraglist[fs[0]])
    if fi["ringcount"]==0:
        mol = Chem.RWMol(Chem.CombineMols(Chem.MolFromSmiles("[W]"),fraglist[fs[0]]))
    else:
        mol = Chem.RWMol(Chem.CombineMols(Chem.MolFromSmiles("[W]"),frag_ringopen(fraglist[fs[0]],fi["valence"])))
        mol.AddBond(0,3,fi["bondtypes"][0])
        mol.AddBond(0,2,fi["bondtypes"][1])
    for idx in fs[1:]:
        if mol!= None:
            if idx>=0:
                mol = add_fs_mol(mol,idx,fraglist,infolist)
            else: #this means a ring closure has to be done
                mol = add_fs_mol_rc(mol,idx,fraglist,infolist)
    return mol
    
def get_last_dummy_count(bestmol):
    """ get dummy count for last in list
    """
    total_dummies=0
    if len(bestmol)>0:
        if bestmol[-1] != None:
            for at in bestmol[-1].GetAtoms():
                if at.GetAtomicNum()==0:
                    total_dummies += 1
    return total_dummies
    
def display_last_molecule(bestmol):
    """ helper function to dislpay last mol in bestmol in verbose mode
    mode.
    """
    mol = Chem.RWMol(bestmol[-1])
    mol.RemoveAtom(0)
    for a in mol.GetAtoms():
        a.SetIsAromatic(False) #to prevent errors for intermediary structs
    try:
        Chem.SanitizeMol(mol)
        display(mol)
    except:
        pass
    return

class Inverter:

    def __init__(self,nbits=2048,radius=3,atoms_dict=ATOMSDICT,add_charges=False):
        self.nbits = nbits
        self.radius = radius
        self.atoms_dict = atoms_dict
        self.fraglist = gen_radius_0_molfrags(self.atoms_dict)
        if add_charges==True:
            self.fraglist+=gen_charged_frags()
   
    def f2bit(self,frag):
        """
        this function returns the bit id that corresponds to a certain
        radius 0 identifier.
        """
        bit=-1
        bi={}
        sanitize_aro(frag)
        try:
            AllChem.GetMorganFingerprintAsBitVect(frag,self.radius,nBits=self.nbits,bitInfo=bi)
        except:
            display(frag)
        for k in bi:
            if (0,0) in bi[k]:
                bit=k
        return bit

    def get_filtered_fraglist(self,ifp,fraglist):
        """ retain only fragments that have onbits that are in the target fp
        """
        onbits = list(ifp.GetOnBits())
        filtered_fraglist = [frag for frag in fraglist if self.f2bit(frag) in onbits]
        filtered_infolist = [inspect_frag(frag) for frag in filtered_fraglist]
        return filtered_fraglist,filtered_infolist

    def mol_score(self,mol,ifp,dummypenalty=1):
        """score the molecule for ranking"""
        ref_onbits = set(ifp.GetOnBits())
        if mol == None:
            score = -100
        else:
            mol_onbits = self.get_non_dummy_onbits(mol)
            score = len(mol_onbits)-dummypenalty*len(mol.GetSubstructMatches(Chem.MolFromSmarts("[#0]")))
            for bit in mol_onbits:
                if bit not in ref_onbits:
                    score=-100
        return score

    def get_non_dummy_onbits(self,mol):
        """ input a smiles that can have dummy atoms and get out all ECFP
        on-bits except the ones that involve dummy atoms.
        """
        bi = {}
        if mol!=None:
            bi1 = {}
            sanitize_aro(mol) #to resolve aromatic non aromatic mixtures and impropper aromatization
            AllChem.GetMorganFingerprintAsBitVect(mol,self.radius,nBits=self.nbits,bitInfo=bi1)
            for key in bi1.keys():
                hasdummies = True
                for frag in bi1[key]:
                    if frag[1]>0: #this is because for some reason at radius 0 it cant find the dummy substructure
                        submol=Chem.PathToSubmol(mol,Chem.FindAtomEnvironmentOfRadiusN(mol,frag[1],frag[0]))
                        a = submol.GetSubstructMatches(Chem.MolFromSmarts("[#0]"))
                        b = submol.GetSubstructMatches(Chem.MolFromSmarts("[W]"))
                        if len(a)+len(b)==0:
                            hasdummies = False
                    else: 
                        if mol.GetAtomWithIdx(frag[0]).GetAtomicNum() not in [0,74]:
                            hasdummies = False
                if hasdummies == False:
                    bi[key]=bi1[key]
        return set(bi.keys())
        
    def getECFP(self,mol):
        """ helper function that gives ECFP with the radius and fp length the
        inverter class was initialized with.
        """
        return AllChem.GetMorganFingerprintAsBitVect(mol,self.radius,nBits=self.nbits)
        
    def search_forward(self,backtrack,idxlist,currentpath,filtered_frags,filtered_infolist,ifp,bestmol,bestscore,verbose,dummy_penalty,reduction_tolerance,steps,found):
        """ operation done when not backtracking: rank all the options and pick
        the best one. then decide whether to backtrack.
        """
        scores = []
        winner=None
        for i,fr in enumerate(idxlist):
            scores.append(self.mol_score(build_mol_from_fs([x[0] for x in currentpath]+[fr],filtered_frags,filtered_infolist),ifp,dummypenalty=dummy_penalty))
        scores_argsort = sorted(range(len(scores)), key=scores.__getitem__) # http://stackoverflow.com/questions/3071415/efficient-method-to-calculate-the-rank-vector-of-a-list-in-python
        idxlist_sorted=[idxlist[i] for i in scores_argsort][::-1]
        idxlist_sorted=idxlist_sorted[:len(idxlist_sorted)-scores.count(-100)]
        currentpath.append(idxlist_sorted)
        if len(idxlist_sorted)>0:
            bestscore_current=max(scores) # find something that is more resistan to ring closures.
            modified_bestscore=bestscore-dummy_penalty*max(0,filtered_infolist[idxlist_sorted[-1]]["valence"]-2)
            if bestscore_current>=modified_bestscore-reduction_tolerance: #clarify why some tolerance for reduced extra bits works !!
                if currentpath[-1][0]>0:
                    bestscore=bestscore_current
                bestmol.append(build_mol_from_fs([x[0] for x in currentpath],filtered_frags,filtered_infolist))
                if bestmol[-1]==None:
                    if verbose==True:
                        print("best mol is unparsable. backtracking.")
                    backtrack=True
                else:
                    if verbose==True:
                        display_last_molecule(bestmol)
                    if len(bestmol[-1].GetSubstructMatches(Chem.MolFromSmarts("[#0]")))==0:
                        if len(self.get_non_dummy_onbits(bestmol[-1]))==len(ifp.GetOnBits()):
                            bestmol[-1].GetAtomWithIdx(0).SetIsAromatic(False)
                            try:
                                winner = Chem.GetMolFrags(bestmol[-1], asMols=True)[1]
                                Chem.SanitizeMol(winner)
                                found=True
                                if verbose==True:
                                    print("did it in {}.".format(steps))
                            except Exception as e:
                                if verbose==True:
                                    print("problem molecule, probably unkekulizable")
                                    print(e)
                                backtrack=True

                        else:
                            if verbose==True:
                                print("terminal path but molecule doesnt match")
                            backtrack=True
            else:
                if verbose == True:
                    print("at this point thte sorted winning ids are {}".format(idxlist_sorted))
                    print("need to backtrack bc the amount of onbits went DOWN! it was {} and now its {}".format(modified_bestscore,bestscore_current))
                backtrack=True
        else:
            if verbose==True:
                print("no matching fragments to add. backtracking.")
            backtrack=True
        return backtrack,currentpath,bestmol,bestscore,found,winner
            
    def retrace(self,backtrack,currentpath,bestmol,found,verbose):
        """ opration done when backtracking. retrace in the tree back to an 
        unexplored part. if there are no more options, halt the search.
        """
        bestmol.append(None)
        while len(currentpath[-1])<=1 and len(currentpath)>1:
            currentpath=currentpath[:-1]
        if len(currentpath)==1:
            if verbose==True:
                print("just had to crawl back all the way to a new radius 0 fragment start.")
        currentpath=currentpath[:-1]+[currentpath[-1][1:]]
        if len(currentpath[0])==0 and found != "failure":
            if verbose==True:
                print("it failed. i am sorry")
            found = "failure"
        else:
            backtrack=False
        return backtrack,currentpath,bestmol,found
    

    def invert(self,ifp,max_steps=5000,verbose=False,max_atoms=50,dummy_penalty=1,reduction_tolerance=10,max_dummies=10):
        """ main inversion routine.
        """
        filtered_frags,filtered_infolist = self.get_filtered_fraglist(ifp,self.fraglist)
        if verbose==True:
            print("these are the filtered frags")
            img = Draw.MolsToGridImage(filtered_frags,molsPerRow=8)
            display(img)
            for i,fra in enumerate(filtered_frags):
                print("frag",i)
                print([b.GetBondType() for b in fra.GetAtomWithIdx(0).GetBonds()])
        idxlist = [x for x in range(len(filtered_frags))]
        #currentpath = [list(idxlist)] # order for valence after this
        currentpath = [[]]
        for curr_valence in range(7):
            currentpath[0] += [idx for idx in idxlist if filtered_infolist[idx]["valence"] == curr_valence]
        idxlist+=[-y-1 for y in range(5)]
        steps = 0
        bestmol = []
        found=False
        backtrack=False
        bestscore=-9999
        winner=None
        while max_steps>steps and found==False:
            total_dummies=get_last_dummy_count(bestmol)
            if backtrack==False and len([x[0] for x in currentpath if x[0]>=0])+1<=max_atoms and len(currentpath[0])>0 and total_dummies<=max_dummies:
                backtrack,currentpath,bestmol,bestscore,found,winner = self.search_forward(backtrack,idxlist,currentpath,filtered_frags,filtered_infolist,ifp,bestmol,bestscore,verbose,dummy_penalty,reduction_tolerance,steps,found)
            else:
                backtrack,currentpath,bestmol,found=self.retrace(backtrack,currentpath,bestmol,found,verbose)
            if verbose == True: 
                print([c[0] for c in currentpath if len(c)>0])
                print(bestscore)
            steps += 1
        return winner, steps

if __name__ == "__main__":
    #display is used for rendering molecules if the function is used in a jupyter notebook.
    #if not the functino is passed. 
    def display(x): 
        pass
    inv1=Inverter()
    inv1.invert(inv1.getECFP(Chem.MolFromSmiles("NOON(COF)F")),verbose=True)
    #import cProfile
    #smi="C1CCOC1c1ccccc1c2ccccc2"
    #cProfile.run("inv1.invert(inv1.getECFP(Chem.MolFromSmiles(smi)),verbose=True)")
