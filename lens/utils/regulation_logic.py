from __future__ import absolute_import, division, print_function

from parsimonious.grammar import Grammar
from parsimonious.nodes import NodeVisitor

grammar = Grammar(
    """
    rule = if one_molecule
    one_molecule = text ws? and? or? not?

    text = ~"[A-Z 0-9]*"i

    if = ws? "IF" ws?
    and = ws? "and" ws?
    or = ws? "or" ws?
    not = ws? "not" ws?
    open  = "("
    close = ")"
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
        pass
    def visit_one_molecule(self, node, visited_children):
        import ipdb; ipdb.set_trace()
        pass
    def visit_text(self, node, visited_children):
        pass
    def visit_if(self, node, visited_children):
        pass
    def visit_and(self, node, visited_children):
        pass
    def visit_or(self, node, visited_children):
        pass
    def visit_not(self, node, visited_children):
        pass
    def visit_open(self, node, visited_children):
        pass
    def visit_close(self, node, visited_children):
        pass
    def visit_ws(self, node, visited_children):
        pass
    def generic_visit(self, node, visited_children):
        # The generic visit method.
        return visited_children or node




# str = "IF (GLCxt or LCTSxt or RIBxt or GLxt or LACxt or PYRxt or SUCCxt or ETHxt or ACxt or FORxt)"
# str2 = "IF not (GLCxt or LCTSxt or RUBxt) and FNR and not GlpR"

str = "IF GLCxt"


str_parsed = grammar.parse(str)
lc = LogicConstructor()
logic = lc.visit(str_parsed)

import ipdb; ipdb.set_trace()
