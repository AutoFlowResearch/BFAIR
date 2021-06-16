import {
  ADD_STUDY,
  REMOVE_STUDY,
  ADD_ASSAY,
  REMOVE_ASSAY,
  SET_ASSAY_TITLE,
  SET_STUDY_TITLE,
} from '../../components/templates/InvestigationFormPage/InvestigationFormPage.constants';
import { Investigation, Study, Assay } from '../model/investigation.model';

const reducer = (
  state: Investigation = new Investigation(),
  action: { type?: string; payload?: any },
) => {
  const { studies } = state;
  const { type, payload } = action;
  switch (type) {
    case ADD_STUDY:
      studies.push(new Study());
      break;
    case REMOVE_STUDY:
      studies.splice(payload, 1);
      break;
    case ADD_ASSAY:
      (payload.study as Study).assays.push(new Assay());
      break;
    case REMOVE_ASSAY:
      const { index } = payload;
      (payload.study as Study).assays.splice(index, 1);
      break;
    case SET_ASSAY_TITLE:
      debugger;
      (payload.assay as Assay).title = payload.title;
      break;
    case SET_STUDY_TITLE:
      (payload.study as Study).title = payload.title;
      break;
  }
  return JSON.parse(JSON.stringify(state));
};

export default reducer;
