import pubchempy
from cobra import Model, Reaction, Metabolite
import cobra
from cobra.io import load_json_model, save_json_model, read_sbml_model, write_sbml_model, load_matlab_model, save_matlab_model

class ModelHandler(object):
    def __init__(self): 
        self.model = Model()
        self.metabolites = {} # A dictionary object that consists of metabolite abbreviation as key which plays the role of metabolite   
        # id as well and a cobra.Metabolite object as a value, this dictionary is filled by the metabolite file provided from the user or 
        # direct assignment with the same <key, value> pairs...            
        self.reactions = {} # A dictionary object that consists of reaction id as key and a cobra.Reaction object as a value,
        # this dictionary is filled by the reactions file provided from the user or direct assignment with the same <key, value> pairs...
        self.special_metabolites = ["NAD", "NADH", "ATP", "AMP", "ADP", "PI"] # Those metabolites is best searched by their common
        #  abbreviation than their names...

    # Preconditions: A path to a file consists of semi-colon separated values without spaces i.e,Abbreviation;metabolite name;compartment
    # the abbreviation must be the same in the included reaction in the reactions file. All the values in this file is important and
    # cannot be left empty...
    # Postconditions: The class variable metabolites will be filled by the given metabolites in the file... 
    def Set_Metabolites_From_File(self, file_path):
        with open(file_path) as f:
            lines = f.readlines()
        
        for line in lines:
            line_elements = line.strip().split(';')
            if line_elements[0].upper() in self.special_metabolites:
                if line_elements[0].upper() == "PI":
                    metabolite_formula = pubchempy.get_compounds("Phosphate", "name")[0].molecular_formula
                else:
                    metabolite_formula = pubchempy.get_compounds(line_elements[0], "name")[0].molecular_formula
            else:
                metabolite_formula = pubchempy.get_compounds(line_elements[1], "name")[0].molecular_formula
            
            metabolite = Metabolite(id=line_elements[0], name=line_elements[1], formula=metabolite_formula, compartment=line_elements[2])
            self.metabolites[line_elements[0]] = metabolite
    

    # Preconditions: The metabolites and reactions objects are filled with data, a provided model name and optionally a model id, if not
    # procided the model name will be concatenated with '_1' to act as the model_id
    # Postconditions: The model object will be filled with reactions to build a cobra.Model object with the required data...
    def Set_Model(self, model_name, model_id=None):
        if model_id == None:
            model_id = model_name + '_1'
        
        self.model = Model(id_or_model=model_id, name=model_name)
        all_reactions = []
        for key in self.reactions.keys():
            all_reactions.append(self.reactions[key])

        self.model.add_reactions(all_reactions)

    def Get_S(self):
        return cobra.util.array.create_stoichiometric_matrix(self.model)
    
    # Preconditions: The model object is feeded by reactions and provided file path, as well to file type if not specified an sbml
    # file will be generated...
    # Postconditios: A file with the specified file type is generated...
    def Save_Model(self, file_path, file_type='sbml'):
        if file_type.lower() == 'json':
            save_json_model(self.model, file_path)
        elif file_type.lower() == 'matlab':
            save_matlab_model(self.model, file_path)
        else:
            write_sbml_model(self.model, file_path)

    #Preconditions: A path to a file consists of semi-colon separated values without spaces i.e,Reaction_id;Reaction_name;the flux;
    # lower_bound;upper_bound make sure that the flux metabolites is consistent with the metabolite abbreviations from the metabolites file.
    # The lower and upper bounds are optional but its better to provide them, the empty bounds will be 0, and 1000 if the fulx is not bidirectional
    # and -1000 and 1000 if  not. The other information is important and cannot be left empty...
    # Postconditions: The class variable reactions will be filled by the given reactions in the file... 
    def Set_Reactions_From_File(self, file_path):
        with open(file_path) as f:
            lines = f.readlines()
        
        for line in lines:
            reactants = []
            products = []
            line_elements = line.strip().split(';')
            flux = line_elements[2]
            reaction = Reaction(id=line_elements[0])
            reaction.name = line_elements[1]

            # Check whether the flux bidirectional or not...
            if flux.find('<->') != -1:
                # Is all the information provided or not and deals with empty bounds...
                if len(line_elements) == 5:
                    if line_elements[3] == '': # No lower bounds...
                        reaction.lower_bound = -1000.0
                    elif line_elements[4] == '': # No upper bounds...
                        reaction.upper_bound = 1000.0
                    else:
                        reaction.lower_bound = float(line_elements[3])
                        reaction.upper_bound = float(line_elements[4])

                elif len(line_elements) == 4: # Only the lower bounds are provided...
                    reaction.lower_bound = float(line_elements[3])
                    reaction.upper_bound = 1000.0
                
                else:
                    reaction.lower_bound = -1000.0
                    reaction.upper_bound = 1000.0
                
                # The flux has reactants and products or not...
                if len(flux.split(' <-> ')) == 2: # if true then the flux has the two components...
                    reactants = flux.split(' <-> ')[0]
                    products = flux.split(' <-> ')[1]
                    reactants = reactants.split(' + ')
                    products = products.split(' + ')

                    # These loops adds one metabolite per turn and they parse each metabolite to get their stoichiometric coeffiecients,
                    # The first loop designed for the reactants and the second for the products...
                    for reactant in reactants:
                        count = 0
                        stoichiometry = ''
                        for char in reactant:
                            if char in ['1','2','3','4','5','6','7','8','9','.']:
                                stoichiometry += char
                                count += 1
                            else:
                                break
                        if stoichiometry == '':
                            stoichiometry = '1'
                        reaction.add_metabolites({
                            self.metabolites[reactant[count:]]: -(float(stoichiometry))
                        })

                    for product in products:
                        count = 0
                        stoichiometry = ''
                        for char in product:
                            if char in ['1','2','3','4','5','6','7','8','9','.']:
                                stoichiometry += char
                                count += 1
                            else:
                                break
                        if stoichiometry == '':
                            stoichiometry = '1'
                        reaction.add_metabolites({
                            self.metabolites[product[count:]]: float(stoichiometry)
                        })

                # This one if the flux has either reactants or products...
                elif len(flux.split(' <-> ')) == 1:
                    if flux.index('<->') == 0: # This means it has products...
                        products = flux.split('<-> ')[1]

                        count = 0
                        stoichiometry = ''
                        for char in products:
                            if char in ['1','2','3','4','5','6','7','8','9','.']:
                                stoichiometry += char
                                count += 1
                            else:
                                break
                        if stoichiometry == '':
                            stoichiometry = '1'
                        reaction.add_metabolites({
                            self.metabolites[products[count:]]: float(stoichiometry)
                        })

                    else:
                        reactants = flux.split(' <->')[0]   

                        count = 0
                        stoichiometry = ''
                        for char in reactants:
                            if char in ['1','2','3','4','5','6','7','8','9','.']:
                                stoichiometry += char
                                count += 1
                            else:
                                break
                        if stoichiometry == '':
                            stoichiometry = '1'
                        reaction.add_metabolites({
                            self.metabolites[reactants[count:]]: -(float(stoichiometry))
                        })
            else: # if the flux is not bidirectional and the other procedures are the same like the first ones...
                if len(line_elements) == 5:
                    if line_elements[3] == '':
                        reaction.lower_bound = 0.0
                    elif line_elements[4] == '':
                        reaction.upper_bound = 1000.0
                    else:
                        reaction.lower_bound = float(line_elements[3])
                        reaction.upper_bound = float(line_elements[4])

                elif len(line_elements) == 4:
                    reaction.lower_bound = float(line_elements[3])
                    reaction.upper_bound = 1000.0
                
                else:
                    reaction.lower_bound = 0.0
                    reaction.upper_bound = 1000.0

                if len(flux.split(' -> ')) == 2:
                    reactants = flux.split(' -> ')[0]
                    products = flux.split(' -> ')[1]
                    reactants = reactants.split(' + ')
                    products = products.split(' + ')

                    for reactant in reactants:
                        count = 0
                        stoichiometry = ''
                        for char in reactant:
                            if char in ['1','2','3','4','5','6','7','8','9','.']:
                                stoichiometry += char
                                count += 1
                            else:
                                break
                        if stoichiometry == '':
                            stoichiometry = '1'
                        reaction.add_metabolites({
                            self.metabolites[reactant[count:]]: -(float(stoichiometry))
                        })

                    for product in products:
                        count = 0
                        stoichiometry = ''
                        for char in product:
                            if char in ['1','2','3','4','5','6','7','8','9','.']:
                                stoichiometry += char
                                count += 1
                            else:
                                break
                        if stoichiometry == '':
                            stoichiometry = '1'
                        reaction.add_metabolites({
                            self.metabolites[product[count:]]: float(stoichiometry)
                        })

                elif len(flux.split(' -> ')) == 1:
                    if flux.index('->') == 0:
                        products = flux.split('-> ')[1]

                        count = 0
                        stoichiometry = ''
                        for char in products:
                            if char in ['1','2','3','4','5','6','7','8','9','.']:
                                stoichiometry += char
                                count += 1
                            else:
                                break
                        if stoichiometry == '':
                            stoichiometry = '1'
                        reaction.add_metabolites({
                            self.metabolites[products[count:]]: float(stoichiometry)
                        })
                    else:
                        reactants = flux.split(' ->')[0]

                        count = 0
                        stoichiometry = ''
                        for char in reactants:
                            if char in ['1','2','3','4','5','6','7','8','9','.']:
                                stoichiometry += char
                                count += 1
                            else:
                                break
                        if stoichiometry == '':
                            stoichiometry = '1'
                        reaction.add_metabolites({
                            self.metabolites[reactants[count:]]: -(float(stoichiometry))
                        })
                
            
            self.reactions[line_elements[0]] = reaction 
