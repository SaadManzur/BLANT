import numpy as np

from matplotlib import pyplot as plt

from sklearn.metrics import auc
from sklearn.metrics import ndcg_score
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.preprocessing import OneHotEncoder
from sklearn.metrics import precision_recall_curve

from utils.compiled_result import CompiledResult

def get_auroc(result, out):
    
    assert isinstance(result, CompiledResult)
    
    y_pred = [1]*len(result)
    
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    
    fprs, tprs, threshold = roc_curve(result.Ts, result.confidences)

    ax.plot(fprs, tprs)
    
    plt.savefig(f'{out}_roc.png')
    
    return roc_auc_score(result.Ts, result.confidences)

def get_aupr(result, out):
    
    assert isinstance(result, CompiledResult)
    
    y_pred = [1]*len(result)
    
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    
    precision, recall, threshold = precision_recall_curve(result.Ts, result.confidences)
    
    ax.plot(recall, precision)
    
    plt.savefig(f'{out}_pr.png')
    
    return auc(recall, precision)

def get_ndcg(result):
    
    assert isinstance(result, CompiledResult)
    
    y_pred = [[1-confidence, confidence] for confidence in result.confidences]
    
    encoder = OneHotEncoder(sparse=False)

    reshaped_ts = np.array(result.Ts).reshape((len(result), 1))
    
    y_true = encoder.fit_transform(reshaped_ts)
    
    return ndcg_score(y_true, y_pred)