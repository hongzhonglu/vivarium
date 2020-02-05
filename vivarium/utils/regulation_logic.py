from __future__ import absolute_import, division, print_function

import pprint

from arpeggio import Optional, ZeroOrMore, OneOrMore, EOF, ParserPython, Kwd, RegExMatch


pp = pprint.PrettyPrinter(indent=4)

def convert_number(s):
    try:
        return float(s)
    except:
        return None


# each function returns a tuple
def symbol(): return RegExMatch(r'[a-zA-Z0-9.\-\_]+')  # TODO -- surplus can be evaluated if there is a threshold, ignored for now
def tuple_key(): return Kwd("("), symbol, Kwd(","), symbol, Kwd(")")
def operator(): return [Kwd('>'), Kwd('<')]
def compare(): return [symbol, tuple_key], Optional(operator, symbol)
def group(): return Kwd("["), logic, Kwd("]")
def term(): return Optional(Kwd("not")), [compare, group]
def logic(): return term, ZeroOrMore([Kwd("and"), Kwd("or")], term)
def rule(): return Kwd("if"), logic, EOF

def evaluate_symbol(tree, state):
    symbol = tree
    value = state.get(symbol.value)
    if value is None:
        value = convert_number(symbol.value)

    return value

def evaluate_tuple_key(tree, state):
    tuple_key = (tree[1].value, tree[3].value)
    value = state.get(tuple_key)

    return value

def evaluate_compare(tree, state):
    if tree[0].rule_name == 'tuple_key':
        first = evaluate_tuple_key(tree[0], state)
    else:
        first = evaluate_symbol(tree[0], state)

    if len(tree) == 1:
        return first
    else:
        operator = tree[1]
        last = evaluate_symbol(tree[2], state)
        if operator.value == '<':
            return first < last
        else:
            return first > last

def evaluate_group(tree, state):
    logic = tree[1]
    return evaluate_logic(logic, state)

def evaluate_term(tree, state):
    invert = False
    value = False

    if tree[0].value == 'not':
        invert = True
        tree = tree[1:]

    if tree[0].rule_name == 'group':
        value = evaluate_group(tree[0], state)
    elif tree[0].rule_name == 'compare':
        value = evaluate_compare(tree[0], state)

    if invert:
        value = not value

    return value

def evaluate_logic(tree, state):
    head = evaluate_term(tree[0], state)

    if len(tree) > 1:
        tail = evaluate_logic(tree[2:], state)
        operation = tree[1].value

        if operation == 'and':
            head = head and tail
        elif operation == 'or':
            head = head or tail

    return head

def evaluate_rule(tree, state):
    return evaluate_logic(tree[1], state)


# make parser based on "rule", defined in grammar above
rule_parser = ParserPython(rule)

def build_rule(expression):
    # type: (str) -> Callable[Dict[str, bool], bool]

    '''
    Accepts a string representing a logical statement about the presence or absence of
    various molecular entities relevant to regulation, and returns a function that
    evaluates that logic with respect to actual values for the various symbols. 
    '''
    tree = rule_parser.parse(expression)
    # print("{}: {}".format(expression, tree))

    def logic(state):
        return evaluate_rule(tree, state)

    return logic

def test_arpeggio():

    # test logic with sets
    test = "if not [GLCxt or LCTSxt or RUBxt] and FNR and not GlpR"
    state_false = {'GLCxt': True, 'LCTSxt': False, 'RUBxt': True, 'FNR': True, 'GlpR': False}
    state_true = {'GLCxt': False, 'LCTSxt': False, 'RUBxt': False, 'FNR': True, 'GlpR': False}
    run_rule = build_rule(test)
    assert run_rule(state_false) == False
    assert run_rule(state_true) == True

    # test regulation logic in set
    test = "if not [FDP > 10 or F6P]"
    state_false = {'FDP': 20, 'F6P': False}
    state_true = {'FDP': 5, 'F6P': False}
    run_rule = build_rule(test)
    assert run_rule(state_false) == False
    assert run_rule(state_true) == True

    # regulation with thresholds
    test = "if not external_glc > 0.1"
    state_false = {'external_glc': 0.2}
    state_true = {'external_glc': 0.01}
    run_rule = build_rule(test)
    assert run_rule(state_false) == False
    assert run_rule(state_true) == True

    # test tuple keys
    test = "if not (external, glc) > 0.1"
    state_false = {('external', 'glc'): 0.2}
    state_true = {('external', 'glc'): 0.01}
    run_rule = build_rule(test)
    assert run_rule(state_false) == False
    assert run_rule(state_true) == True

    print('test passed!')

    return run_rule


if __name__ == '__main__':
    test_arpeggio()
