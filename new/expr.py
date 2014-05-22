#!/usr/bin/env python

import os
import sys
from pyparsing import Regex, Word, Group, operatorPrecedence
from pyparsing import alphas, alphanums, opAssoc

class ExpressionParser:
  def __init__(self):
    self.eq_op="e"
    self.ne_op="n"
    self.like_op="l"
    self.not_like_op="nl"
    self.and_op = "a"
    self.or_op = "o"
    self.in_op = "i"
    self.tags = set()
    self.operators = list()
    self.relational_exps = list()

    def validate_tag(tokens):
      self.tags.update(tokens)
      return tokens[0]

    eq_ops = Regex('==?|is(?i)').setParseAction(lambda t:self.eq_op)
    ne_ops = Regex('!=|is\s+not(?i)').setParseAction(lambda t:self.ne_op)
    like_ops = Regex('like(?i)').setParseAction(lambda t:self.like_op)
    not_like_ops = Regex('not\s+like(?i)').setParseAction(lambda t:self.not_like_op)
    and_ops = Regex('&&|and(?i)').setParseAction(lambda t:self.and_op)
    or_ops = Regex('\|\||or(?i)').setParseAction(lambda t:self.or_op)

    #the order of "ne_ops" and "eq_ops" matters
    relational_ops = ne_ops | eq_ops | like_ops | not_like_ops
    logical_ops = and_ops | or_ops

    tag = Word(alphas, alphanums + '_').setParseAction(validate_tag)
    number = Regex('-?[0-9]+(.[0-9]+)?').setParseAction(lambda t:'"' + t[0] + '"')
    literal = Regex('".*"').setParseAction(lambda t:t[0][1:-1])
    null_value = Regex('null(?i)').setParseAction(lambda v:"NULL")
    term = null_value | tag | number | literal

    relational_exp = Group(term + relational_ops + term)
    self.logical_exp = operatorPrecedence(relational_exp, [(logical_ops, 2, opAssoc.LEFT, collect_data)], lpar=("\0"), rpar="\0")

  

    self.logical_exp.parseString('ABC=def || DGc==324.3 AND bcd != -34 OR "sdfdf" is te && FDF is   not null', parseAll=True)
    print result
    
    print self.logical_exp.parseString('ABC=def', parseAll=True)
    print self.logical_exp.parseString('ADB=def is null', parseAll=True)
  
  def collect_data
  def simple_optimize(exp):
    return 1

def main():
  ep = ExpressionParser()

if __name__ == "__main__":
  main()
