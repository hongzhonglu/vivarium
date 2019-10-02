from __future__ import absolute_import, division, print_function

import pprint

from arpeggio import Optional, ZeroOrMore, OneOrMore, EOF, ParserPython, Kwd, RegExMatch


pp = pprint.PrettyPrinter(indent=4)

ignored_statements = ['action is complex', '']

def symbol(): return Optional(Kwd("surplus")), RegExMatch(r'[a-zA-Z0-9\[\]\-\_]+')  # TODO -- surplus can be evaluated if there is a threshold, ignored for now
def group(): return Kwd("("), logic, Kwd(")")
def term(): return Optional(Kwd("not")), [symbol, group]
def logic(): return term, ZeroOrMore([Kwd("and"), Kwd("or")], term)
def rule(): return Optional(Kwd("active")), Kwd("IF"), logic, EOF

def evaluate_symbol(tree, env):
    symbol = tree[0]
    if symbol == 'surplus':
        symbol = tree[1]
    value = env.get(symbol.value)

    return value

def evaluate_group(tree, env):
    logic = tree[1]
    return evaluate_logic(logic, env)

def evaluate_term(tree, env):
    invert = False
    value = False

    if tree[0].value == 'not':
        invert = True
        tree = tree[1:]

    if tree[0].rule_name == 'group':
        value = evaluate_group(tree[0], env)
    elif tree[0].rule_name == 'symbol':
        value = evaluate_symbol(tree[0], env)

    if invert:
        value = not value

    return value

def evaluate_logic(tree, env):
    head = evaluate_term(tree[0], env)

    if len(tree) > 1:
        tail = evaluate_logic(tree[2:], env)
        operation = tree[1].value

        if operation == 'and':
            head = head and tail
        elif operation == 'or':
            head = head or tail

    return head

def evaluate_rule(tree, env):
    if tree[0].value == 'active':
        tree = tree[1:]
    return evaluate_logic(tree[1], env)


rule_parser = ParserPython(rule)

def null_fun(env):
    return None

def build_rule(expression):
    # type: (str) -> Callable[Dict[str, bool], bool]

    '''
    Accepts a string representing a logical statement about the presence or absence of
    various molecular entities relevant to regulation, and returns a function that
    evaluates that logic with respect to actual values for the various symbols. 
    '''
    if expression in ignored_statements:
        return null_fun
    else:
        tree = rule_parser.parse(expression)
        def logic(env):
            return evaluate_rule(tree, env)

        return logic

def test_arpeggio():
    test = "IF not (GLCxt or LCTSxt or RUBxt) and FNR and not GlpR"
    state_false = {'GLCxt': True, 'LCTSxt': False, 'RUBxt': True, 'FNR': True, 'GlpR': False}
    state_true = {'GLCxt': False, 'LCTSxt': False, 'RUBxt': False, 'FNR': True, 'GlpR': False}
    run_rule = build_rule(test)
    assert run_rule(state_false) == False
    assert run_rule(state_true) == True

    # test surplus in statement
    test = "active IF not (surplus FDP or F6P)"
    state_false = {'FDP': True, 'F6P': False}
    state_true = {'FDP': False, 'F6P': False}
    run_rule = build_rule(test)
    assert run_rule(state_false) == False
    assert run_rule(state_true) == True

    test = 'action is complex'
    run_rule = build_rule(test)
    assert run_rule(state_false) == None

    return run_rule


if __name__ == '__main__':
    test_arpeggio()
