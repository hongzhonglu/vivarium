from __future__ import absolute_import, division, print_function

from parsimonious.grammar import Grammar
from parsimonious.nodes import NodeVisitor

grammar = Grammar(
    """
    rule = if set_mols*
    set_mols = open? one_molecule* close?
    one_molecule = operation? text
    text = ~"[A-Za-z0-9]*"i
    if = "IF" ws
    open  = "("
    close = ")"
    operation = or? and? not?
    or = ws "or" ws?
    and = ws "and" ws?
    not = ws "not" ws?
    ws = ~"\s*"
    """)

class RegulatoryLogic(object):
    def __init__(self):
        self.logic_constructor = LogicConstructor()

    def get_logic(self, logic_str):
        '''
        Make a logic function from a string
        Args:
            logic_str (str)
        '''
        logic_parsed = grammar.parse(logic_str)
        logic_function = self.logic_constructor.visit(logic_parsed)
        return logic_function

class LogicConstructor(NodeVisitor):
    '''
    Make a logic function from a parsed expression.
    Args:
        - node: The node we're visiting
        - visited_children: The results of visiting the children of that node, in a list
    '''
    def visit_rule(self, node, visited_children):
        if_statement, sets_mols, = visited_children
        rule_string = 'if '
        for set_mols in sets_mols:
            for mol_operation in set_mols:
                operation, mol = mol_operation
                if operation:
                    mol_dict = '{} dict.get({}, False)'.format(operation, mol)
                else:
                    mol_dict = 'dict.get({}, False)'.format(operation, mol)
                rule_string = rule_string + mol_dict + ' '
        return rule_string
    def visit_set_mols(self, node, visited_children):
        open, molecules, close = visited_children

        # import ipdb;
        # ipdb.set_trace()


        return molecules
    def visit_one_molecule(self, node, visited_children):
        operation_list, mol_id = visited_children
        operation = operation_list[0]  # TODO -- can this be take out of list before?
        return (operation, mol_id)
    def visit_text(self, node, visited_children):
        return (node.text)
    def visit_operation(self, node, visited_children):
        or_oper, and_oper, not_oper = visited_children
        if isinstance(or_oper, list):
            return ('or')
        elif isinstance(and_oper, list):
            return ('and')
        elif isinstance(not_oper, list):
            return ('not')
    def visit_or(self, node, visited_children):
        return True
    def visit_and(self, node, visited_children):
        return True
    def visit_not(self, node, visited_children):
        return True
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




# str = "IF (GLCxt or LCTSxt or RIBxt or GLxt or LACxt or PYRxt or SUCCxt or ETHxt or ACxt or FORxt)"
str = "IF not (GLCxt or LCTSxt or RUBxt) and FNR and not GlpR"

# str = "IF GLCxt or LCTSxt"

str_parsed = grammar.parse(str)
lc = LogicConstructor()
logic = lc.visit(str_parsed)

# import ipdb; ipdb.set_trace()



def test_rule(dict):
    return (dict.get('GLCxt', False) or dict.get('LCTSxt', False)) and dict.get('FNR', False) and not dict.get('GlpR', False)


state = {
    'GLCxt': False,
    'LCTSxt': True,
    'FNR': True,
    'GlpR': True,
}

result = test_rule(state)

import ipdb; ipdb.set_trace()
