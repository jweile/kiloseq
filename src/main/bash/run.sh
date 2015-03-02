#!/bin/bash
nohup lib/master.R $@ >kiloseq.out 2>kiloseq.err &
