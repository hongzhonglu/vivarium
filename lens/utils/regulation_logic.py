from __future__ import absolute_import, division, print_function

import pprint
pp = pprint.PrettyPrinter(indent=4)

from parsimonious.grammar import Grammar
from parsimonious.nodes import NodeVisitor


nodes = Grammar(
	"""
	
	""")




# nodes = Grammar(
#     """
#     rule = active? if logic (operator logic)*
#     logic = group / operation
# 	group = not? "(" surplus? operation ")"
#     operation = node (operator logic)*
#     node = not? surplus? text
#     operator = or / and

#     or = ws "or" ws
#     and = ws "and" ws
#     not = "not" ws

#     text = ~"[A-Za-z0-9-\[\]]+"i

#     if = "IF" ws
#     active = "active" ws
#     surplus = "surplus" ws
#     open  = "("
#     close = ")"
#     ws = ~"\s+"
#     """)


simplify = Grammar(
    """
    rule = active? if set_mols+
    set_mols = not? open? one_molecule+ close?
    one_molecule = surplus? operation? not? text
    text = ~"[A-Za-z0-9-\[\]]+"i
    if = "IF" ws
    active = "active" ws
    surplus = "surplus" ws
    open  = "("
    close = ")"
    operation = or / and
    or = ws "or" ws
    and = ws "and" ws
    not = "not" ws
    ws = ~"\s+"
    """)


class Key(object):
    def __init__(self, key):
        self.key = key

    def __repr__(self):
        return "Key({})".format(self.key)


class GenericConstructor(NodeVisitor):
    def generic_visit(self, node, visited_children):
        if node.expr_name:
            value = [Key(node.expr_name), node.text]
            for child in visited_children:
                if isinstance(child, list):
                    if child:
                        if isinstance(child[0], Key):
                            value.append(child)
                        else:
                            for grandchild in child:
                                value.append(grandchild)
                else:
                    value.append(child)
            return value
        else:
            return visited_children


def pull_values(tree):
    key = tree[0]
    string = tree[1]
    children = tree[2:] if len(tree) > 2 else []
    
    return key, string, children

def evaluate_logic(outcome, operation, inverse, value):
    value = not value if inverse else value

    if operation:
        if operation == 'and':
            outcome = outcome and value
        elif operation == 'or':
            outcome = outcome or value
    else:
        outcome = value

    return outcome

def evaluate_mol(tree, env):
    key, string, children = pull_values(tree)

    operation = None
    inverse = None
    if children[0][0].key == 'surplus':
        children = children[1:]
    if children[0][0].key == 'operation':
        operation = children[0][2][0].key
        children = children[1:]
    if children[0][0].key == 'not':
        inverse = True
        children = children[1:]
    value = env.get(children[0][1])

    print('key({}): {} - {}'.format(key, string, children))
    print('{}operation({}): {}'.format('not ' if inverse else '', operation, value))

    return operation, inverse, value

def evaluate_set_mol(tree, env):
    key, string, children = pull_values(tree)
    
    inverse = children[0][0].key == 'not'
    mols = [child for child in children if child[0].key == 'one_molecule']
    outcome = True
    
    first_operation = None
    first_inverse = None
    for index, mol in enumerate(mols):
        operation, inverse, value = evaluate_mol(mol, env)
        if index == 0:
            first_operation = operation
            first_inverse = inverse
        outcome = evaluate_logic(outcome, operation, inverse, value)
        print('outcome: {}'.format(outcome))

    print('key({}): {} - {}'.format(key, string, children))
    print('{}operation({}): {}'.format('not ' if first_inverse else '', first_operation, value))

    return operation, inverse, outcome

def evaluate_rule(tree, env):
    key, string, children = pull_values(tree)

    # ingore active and if keys for now
    set_mols = [child for child in children if child[0].key == 'set_mols']
    outcome = True

    for set_mol in set_mols:
        operation, inverse, value = evaluate_set_mol(set_mol, env)
        outcome = evaluate_logic(outcome, operation, inverse, value)

    return outcome
    




grammar = Grammar(
    """
    rule = active? if set_mols+
    set_mols = operation? open? one_molecule+ close?
    one_molecule = surplus? operation? operation? text
    text = ~"[A-Za-z0-9-\[\]]*"i
    if = "IF" ws
    active = "active" ws
    surplus = "surplus" ws
    open  = "("
    close = ")"
    operation = or / and / not
    or = ws "or" ws?
    and = ws "and" ws?
    not = ws "not" ws?
    ws = ~"\s*"
    """)

class RegulatoryLogic(object):
    def __init__(self):
        self.logic_constructor = LogicConstructor()

    def get_logic_function(self, logic_str):
        '''
        Make a logic function from a string
        Args:
            logic_str (str)
        Returns: logic_function (function) that takes in a dict with boolean values {mol_id: bool},
            and evaluates it according to the parsed expression

        '''

        try:
            logic_parsed = grammar.parse(logic_str)
            return self.logic_constructor.visit(logic_parsed)
        except:
            def fun(dict):
                return None
            return fun

class LogicConstructor(NodeVisitor):
    '''s
    Make a logic function from a parsed expression.
    Args:
        - node: The node we're visiting
        - visited_children: The results of visiting the children of that node, in a list
    Returns:
        - a logic function that takes in a dict with boolean values {mol_id: bool},
            and evaluates it according to the parsed expression
    '''
    def visit_rule(self, node, visited_children):
        active_statement, if_statement, sets_mols, = visited_children
        rule_string = ''
        for logic_set in sets_mols:
            set_operation = logic_set[0][0]
            set_mols = logic_set[1]
            in_set = logic_set[2]

            rule_string = rule_string + ' ' + set_operation + ' '
            if in_set:
                rule_string = rule_string + '('
            for mol_operation in set_mols:
                operations, mol = mol_operation
                operation1, operation2 = operations

                if operation1 and operation2:
                    mol_dict = "{} {} dict.get('{}', False)".format(operation1, operation2, mol)
                elif operation1:
                    mol_dict = "{} dict.get('{}', False)".format(operation1, mol)
                else:
                    mol_dict = "dict.get('{}', False)".format(mol)
                rule_string = rule_string + mol_dict + ' '
            rule_string = rule_string[:-1]
            if in_set:
                rule_string = rule_string + ')'

        def logic_function(dict):
            return eval(rule_string)

        return logic_function

    def visit_set_mols(self, node, visited_children):
        operation, open_set, molecules, close_set = visited_children
        in_set = False
        if not isinstance(operation, list):
            operation =['']
        if isinstance(open_set, list) or isinstance(close_set, list):
            assert (isinstance(open_set, list) and isinstance(close_set, list))
            in_set = True
        return (operation, molecules, in_set)

    def visit_one_molecule(self, node, visited_children):
        surplus, operation_list1, operation_list2, mol_id = visited_children

        operation1 = ''
        operation2 = ''
        if isinstance(operation_list1, list):
            operation1 = operation_list1[0]
        if isinstance(operation_list2, list):
            operation2 = operation_list2[0]
        return ([operation1, operation2], mol_id)

    def visit_operation(self, node, visited_children):
        ''' if there is an operation, pull the string out of the list and return it'''
        oper_list = visited_children
        if isinstance(oper_list, list):
            return oper_list[0]
    def visit_text(self, node, visited_children):
        return (node.text)
    def visit_or(self, node, visited_children):
        return ('or')
    def visit_and(self, node, visited_children):
        return ('and')
    def visit_not(self, node, visited_children):
        return ('not')
    def visit_surplus(self, node, visited_children):  # TODO -- base whether surplus on a threshold?
        pass
    def visit_active(self, node, visited_children):
        pass
    def visit_if(self, node, visited_children):
        pass
    def visit_open(self, node, visited_children):
        return ('(')
    def visit_close(self, node, visited_children):
        return (')')
    def visit_ws(self, node, visited_children):
        pass
    def generic_visit(self, node, visited_children):
        # The generic visit method.
        return visited_children or node



# # str = "IF GLCxt or LCTSxt"
# # str = "IF (GLCxt or LCTSxt or RIBxt or GLxt or LACxt or PYRxt or SUCCxt or ETHxt or ACxt or FORxt)"
# # str = "IF not (GLCxt or LCTSxt or RUBxt)"
# # str = "IF not (GLCxt or LCTSxt or RUBxt) and FNR and not GlpR"
# # str = "IF not (GLCxt or LCTSxt or RUBxt) and FNR and GlpR"
# # str = "active IF not (OXYGEN-MOLECULE[e])"
# # str = "active IF not (surplus FDP or F6P)"
# str = 'action is complex'
#
#
# rc = RegulatoryLogic()
# logic_function = rc.get_logic_function(str)
#
# state = {
#     'GLCxt': True,
#     # 'LCTSxt': True,
#     # 'RUBxt': False,
#     # 'FNR': True,
#     # 'GlpR': False,
#     'OXYGEN-MOLECULE[e]': False,
#     # 'FDP': False,
#     # 'F6P': False,
# }
# #
# result = logic_function(state)
# print("RESULT: {}".format(result))


def test_parsing():
    # test = "IF not (GLCxt or LCTSxt or RUBxt) and FNR and not GlpR"
    test = "IF not (GLCxt or LCTSxt or RUBxt) and FNR and not GlpR"
    state_false = {'GLCxt': True, 'LCTSxt': False, 'RUBxt': True, 'FNR': True, 'GlpR': False}
    state_true = {'GLCxt': False, 'LCTSxt': False, 'RUBxt': False, 'FNR': True, 'GlpR': False}

    print(test)

    # tree = nodes.parse(test)
    raw = simplify.parse(test)
    tree = GenericConstructor().visit(raw)

    false = evaluate_rule(tree, state_false)
    print('false: {}'.format(false))
    print('')

    true = evaluate_rule(tree, state_true)

    print('true: {}'.format(true))
    print('')

    return tree
    # return tree, TreeConstructor().visit(tree)

if __name__ == '__main__':
    tree = test_parsing()
    pp.pprint(tree)
