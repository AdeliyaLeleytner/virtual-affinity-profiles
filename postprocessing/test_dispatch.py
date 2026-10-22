"""Portable config/hash dispatch checks with no scientific execution."""
from pathlib import Path
import hashlib,json,tempfile,unittest
from unittest.mock import patch
import dispatch


class DispatchTests(unittest.TestCase):
    def setup_snapshot(self,root,status='complete'):
        inputs=root/'inputs';inputs.mkdir();(inputs/'input.csv').write_text('x\n1\n')
        result=root/'saved';result.mkdir()
        (result/'run_receipt.json').write_text(json.dumps(dict(status=status)))
        (result/'effective_config.json').write_text(json.dumps(dict(datasets=[dict(dataset='toy')])))
        folder=result/'toy__toy__native_genes__all_records';split=folder/'core'/'regime'/'split';split.mkdir(parents=True);(split/'complete.json').write_text('{}')
        files={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for d in [inputs,result] for p in d.rglob('*') if p.is_file()}
        manifest=root/'manifest.json';manifest.write_text(json.dumps(dict(files=files)))
        cfg=root/'config.json';cfg.write_text(json.dumps(dict(postprocessing=dict(input_dir=str(inputs),results=[str(result)],output_dir=str(root/'output'),source_manifest=str(manifest),resamples=0))))
        return cfg,inputs,result

    def test_dispatch_verifies_then_calls_prediction_only_pipeline(self):
        with tempfile.TemporaryDirectory(dir=dispatch.HERE) as temp:
            root=Path(temp);cfg,_,result=self.setup_snapshot(root)
            def fake_process(source,folder,args):
                dest=args.output_dir/source.name/folder.name;dest.mkdir(parents=True)
                (dest/'dataset_state.json').write_text(json.dumps(dict(dataset_complete=True)))
            old_summary=dispatch.postprocess.summary
            try:
                with patch.object(dispatch.postprocess,'process_dataset',side_effect=fake_process) as call,patch.object(dispatch,'support_main') as sup,patch.object(dispatch,'association_main') as asc:
                    dispatch.main(cfg)
                    self.assertEqual(call.call_count,1);self.assertEqual(sup.call_count,1);self.assertEqual(asc.call_count,1)
            finally:dispatch.postprocess.summary=old_summary
            receipt=json.loads((root/'output/dispatch_receipt.json').read_text())
            self.assertTrue(receipt['source_data_unchanged_after_replay']);self.assertFalse(receipt['new_model_training'])

    def test_mismatched_snapshot_hash_rejects_before_any_replay(self):
        with tempfile.TemporaryDirectory(dir=dispatch.HERE) as temp:
            root=Path(temp);cfg,inputs,_=self.setup_snapshot(root)
            (inputs/'input.csv').write_text('changed\n')
            with patch.object(dispatch.postprocess,'process_dataset') as call:
                with self.assertRaisesRegex(RuntimeError,'source hash mismatch'):dispatch.main(cfg)
                call.assert_not_called()

    def test_new_unhashed_source_file_rejects_before_replay(self):
        with tempfile.TemporaryDirectory(dir=dispatch.HERE) as temp:
            cfg,inputs,_=self.setup_snapshot(Path(temp))
            (inputs/"unexpected.csv").write_text("new\n")
            with self.assertRaisesRegex(RuntimeError,"file membership changed"):dispatch.main(cfg)

    def test_incomplete_run_is_not_a_complete_replay(self):
        with tempfile.TemporaryDirectory(dir=dispatch.HERE) as temp:
            cfg,_,_=self.setup_snapshot(Path(temp),status='running')
            with self.assertRaisesRegex(RuntimeError,'Complete source run required'):dispatch.main(cfg)


if __name__=='__main__':unittest.main()
