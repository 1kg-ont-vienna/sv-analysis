import numpy as np

class Family:

    def __init__(self, member) -> None:
        self.members = set([member])
        pass

    def add(self, member):
        self.members.add(member)

    def name(self):
        return '_'.join(sorted(self.members))
    
    def __repr__(self) -> str:
        return self.name()

def merge_family(set1, set2):
    for member in set1:
        set2.add(member)
    return set2

ped_file = np.loadtxt('resources/pedigree.ped', delimiter='\t',dtype=str)[:,1:4]
ped_dict = {line[0]: None for line in ped_file}
for line in ped_file:
    if (line[1] == '0') and (line[2] == '0'):
        if (ped_dict[line[0]] is None):
            ped_dict[line[0]] = Family(line[0])
            continue
        else:
            continue
    
    if (line[1] == '0'):
        if (ped_dict[line[0]] is None):
            if (ped_dict[line[2]] is None):
                ped_dict[line[0]] = Family(line[0])
                ped_dict[line[0]].add(line[2])
                ped_dict[line[2]] = ped_dict[line[0]]
                continue
            else:
                ped_dict[line[2]].add(line[0])
                ped_dict[line[0]] = ped_dict[line[2]]
                continue
        else:
            if (ped_dict[line[2]] is None):
                ped_dict[line[0]].add(line[2])
                ped_dict[line[2]] = ped_dict[line[0]]
                continue
            else:
                members = merge_family(ped_dict[line[0]].members, ped_dict[line[2].members])
                ped_dict[line[0]].members = members
                ped_dict[line[2]] = ped_dict[line[0]]
                continue
    
    if (line[2] == '0'):
        if (ped_dict[line[0]] is None):
            if (ped_dict[line[1]] is None):
                ped_dict[line[0]] = Family(line[0])
                ped_dict[line[0]].add(line[1])
                ped_dict[line[1]] = ped_dict[line[0]]
                continue
            else:
                ped_dict[line[1]].add(line[0])
                ped_dict[line[0]] = ped_dict[line[1]]
                continue
        else:
            if (ped_dict[line[1]] is None):
                ped_dict[line[0]].add(line[1])
                ped_dict[line[1]] = ped_dict[line[0]]
                continue
            else:
                members = merge_family(ped_dict[line[0]].members, ped_dict[line[1].members])
                ped_dict[line[0]].members = members
                ped_dict[line[1]] = ped_dict[line[0]]
                continue
    
    assert (line[1] != '0') and (line[2] != '0')

    try:
        members1 = ped_dict[line[0]].members
    except:
        members1 = set()
    try:
        members2 = ped_dict[line[1]].members
    except:
        members2 = set()
    try:
        members3 = ped_dict[line[2]].members
    except:
        members3 = set()
    
    members = set([line[0], line[1], line[2]])
    members = merge_family(members, members1)
    members = merge_family(members, members2)
    members = merge_family(members, members3)
    ped_dict[line[0]] = Family(line[0])
    ped_dict[line[0]].members = members
    ped_dict[line[1]] = ped_dict[line[0]]
    ped_dict[line[2]] = ped_dict[line[0]]

family_set = set()
for sample in ped_dict:
    family_set.add(ped_dict[sample])
