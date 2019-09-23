from __future__ import absolute_import, division, print_function

from parsimonious.grammar import Grammar
from parsimonious.nodes import NodeVisitor

grammar = Grammar(
    """
    rule = active? if set_mols*
    set_mols = operation? open? one_molecule* close?
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



# str = "IF GLCxt or LCTSxt"
str = "IF (GLCxt or LCTSxt or RIBxt or GLxt or LACxt or PYRxt or SUCCxt or ETHxt or ACxt or FORxt)"
# str = "IF not (GLCxt or LCTSxt or RUBxt)"
# str = "IF not (GLCxt or LCTSxt or RUBxt) and FNR and not GlpR"
# str = "IF not (GLCxt or LCTSxt or RUBxt) and FNR and GlpR"
# str = "active IF not (OXYGEN-MOLECULE[e])"
# str = "active IF not (surplus FDP or F6P)"
# str = 'action is complex'


rc = RegulatoryLogic()
logic_function = rc.get_logic_function(str)

state = {
    'GLCxt': True,
    # 'LCTSxt': True,
    # 'RUBxt': False,
    # 'FNR': True,
    # 'GlpR': False,
    'OXYGEN-MOLECULE[e]': False,
    # 'FDP': False,
    # 'F6P': False,
}
#
result = logic_function(state)
print("RESULT: {}".format(result))
