import { connect } from 'react-redux';
import {
  CONTACTS,
  DESIGN_TYPES,
  FACTOR_TYPES,
  PROTOCOLS,
  PUBLICATIONS,
  STUDY_TYPES,
} from '../../../../../shared/constants';
import { Study } from '../../../../../shared/model/investigation.model';
import { setStudyTitle } from '../../../../../shared/redux/actions';
import Badge from '../../../../atoms/Badge';
import ButtonToggle from '../../../../atoms/ButtonToggle';
import Collapsible from '../../../../atoms/Collapsible';
import DeleteIcon from '../../../../atoms/DeleteIcon';
import Input from '../../../../atoms/Input';
import SelectBar from '../../../../atoms/SelectBar';
import TextArea from '../../../../atoms/TextArea';
import InvestigationFormAssays from '../InvestigationFormAssays';
import './InvestigationFormStudy.style.scss';

const InvestigationStudyTitle = (props: { study: Study; deleteStudy }) => {
  const { study, deleteStudy } = props;
  return (
    <div className='study__title'>
      {study.title || 'New Study'}
      <div className='study__delete'>
        <DeleteIcon onClick={deleteStudy} />
      </div>
      <div className='study__assay-count'>
        <Badge isActive={true}>{study.assays.length} Assay</Badge>
      </div>
    </div>
  );
};

const InvestigationFormStudy = (props: {
  study: Study;
  deleteStudy;
  setStudyTitle?: any;
}) => {
  const { study, deleteStudy, setStudyTitle } = props;
  const handleStudyTitleChange = (title) => {
    setStudyTitle(study, title);
  };
  return (
    <div className='investigation-form__study'>
      <Collapsible
        trigger={
          <InvestigationStudyTitle study={study} deleteStudy={deleteStudy} />
        }
      >
        <Input label='Title' onChange={handleStudyTitleChange} />
        <TextArea label='Description' />
        <SelectBar
          label={'Publications'}
          placeHolderText={'Publication'}
          options={PUBLICATIONS}
          multiple
        />
        <SelectBar
          label={'Contacts'}
          placeHolderText={'Contact'}
          options={CONTACTS}
          multiple
        />
        <ButtonToggle
          label='Study Type'
          items={STUDY_TYPES}
          value={STUDY_TYPES[0]}
          labelKey='studyName'
        />
        <SelectBar
          label={'Design Type'}
          placeHolderText={'Design Type'}
          options={DESIGN_TYPES}
        />
        <Input label='Factor Name' />
        <SelectBar
          label={'Factor Type'}
          placeHolderText={'Factor Type'}
          options={FACTOR_TYPES}
        />
        <SelectBar
          label={'Protocols'}
          placeHolderText={'Protocols'}
          options={PROTOCOLS}
          multiple
        />
        <InvestigationFormAssays study={study} />
      </Collapsible>
    </div>
  );
};

const mapDispatchToProps = (dispatch) => ({
  setStudyTitle: (study, title) => dispatch(setStudyTitle(study, title)),
});

export default connect(null, mapDispatchToProps)(InvestigationFormStudy);
