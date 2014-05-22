#!/usr/bin/env python

import os
import sys
import copy
from pyparsing import Regex, Word, Group, operatorPrecedence
from pyparsing import alphas, alphanums, opAssoc, dblQuotedString

class ExpressionParser:
  def __init__(self):
    self.eq_op="e"
    self.ne_op="n"
    self.like_op="l"
    self.not_like_op="nl"
    self.and_op = "a"
    self.or_op = "o"
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
    number = Regex('-?[0-9]+(.[0-9]+)?').setParseAction(lambda t:[['"' + x + '"' for x in t]])
    literal = copy.copy(dblQuotedString).setParseAction(lambda t:[[ '"' + x[1:-1] + '"' for x in t]])
    null_value = Regex('null(?i)').setParseAction(lambda v:"NULL")
    term = null_value | tag | number | literal

    relational_exp = Group(term + relational_ops + term)
    self.logical_exp = operatorPrecedence(relational_exp, [(logical_ops, 2, opAssoc.LEFT)], lpar=("\0"), rpar="\0")

  

    result = self.logical_exp.parseString('ABC=def || DGc==324.3 AND bcd != -34 OR "sdfdf" is te && FDF is   not null', parseAll=True)
    
    def simple_optimize(exp):
      for i in range(len(exp)):
	  if isinstance(exp[i], list):
	    exp[i] = simple_optimize(exp[i])
      optimized_exp = []
      optimized_exp.append(exp[0])
      
      for i in range(1, len(exp), 2):
	if ((exp[i] == self.or_op and optimized_exp[-1][0] == exp[i + 1][0] and optimized_exp[-1][1] == self.eq_op and exp[i + 1][1] == self.eq_op) 
          or (exp[i] == self.and_op and optimized_exp[-1][0] == exp[i + 1][0] and optimized_exp[-1][1] == self.ne_op and exp[i + 1][1] == self.ne_op)):
	  optimized_exp[-1][2].extend(exp[i + 1][2])
	elif (exp[i] in [self.eq_op, self.ne_op] and not isinstance(exp[i + 1], list) 
	  and (isinstance(optimized_exp[-1], list) or optimized_exp[-1] < exp[i + 1])):
	  last_operand = optimized_exp[-1]
	  optimized_exp[-1] = exp[i + 1]
	  optimized_exp.append(exp[i])
	  optimized_exp.append(last_operand)
	else:	  
	  optimized_exp.append(exp[i])
	  optimized_exp.append(exp[i + 1])
      return optimized_exp

    print simple_optimize(result)

    print simple_optimize(self.logical_exp.parseString('TimeInForce == 6 || 1 == TimeInForce And Test != "abc" And "ghf" is   not Test', parseAll=True).asList())

    print self.logical_exp.parseString('ABC=def', parseAll=True)
    #print self.logical_exp.parseString('ADB=def is null', parseAll=True)
   
def main():
  ep = ExpressionParser()

if __name__ == "__main__":
  main()
