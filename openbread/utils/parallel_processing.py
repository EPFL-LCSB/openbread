# -*- coding: utf-8 -*-
"""
.. module:: openbread
   :platform: Unix, Windows
   :synopsis: Python Multiscale Enzyme Simulator

.. moduleauthor:: openbread team

[---------]

Copyright 2017 Laboratory of Computational Systems Biotechnology (LCSB),
Ecole Polytechnique Federale de Lausanne (EPFL), Switzerland

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

"""


from dill import dumps,loads

def run_dill_encoded(payload):
    fun, args = loads(payload)
    return fun(*args)

def apply_async(pool, fun, args, callback = None):
    payload = dumps((fun, args))
    if callback is None:
        return pool.apply_async(run_dill_encoded, (payload,))
    else:
        return pool.apply_async(run_dill_encoded, (payload,), callback = callback)
