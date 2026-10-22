"""Use tested batched statistics with the unchanged replay/orchestration code."""
import postprocess
from fast_resampling import summary


if __name__=='__main__':
    postprocess.summary=summary
    postprocess.main()
