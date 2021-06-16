import React from 'react';
import {
  CONTACTS,
  ISA_DOCUMENT_LICENSES,
  PUBLICATIONS,
} from '../../../../../shared/constants';
import Collapsible from '../../../../atoms/Collapsible';
import DatePickerInput from '../../../../atoms/DatePickerInput';
import Input from '../../../../atoms/Input';
import SelectBar from '../../../../atoms/SelectBar';
import TextArea from '../../../../atoms/TextArea';

const InvestigationFormGeneral = (props) => {
  return (
    <Collapsible trigger="General Information">
      <Input label="Title" />
      <TextArea label="Description" />
      <div className="row">
        <div className="col-sm-12 col-lg-6 col-xl-4">
          <DatePickerInput label="Submission Date" />
        </div>
        <div className="col-sm-12 col-lg-6 col-xl-4">
          <DatePickerInput label="Public Release Date" />
        </div>
      </div>
      <SelectBar
        label={'Publications'}
        placeHolderText={'Publication'}
        options={PUBLICATIONS}
        multiple
      />
      <SelectBar
        label={'ISA Document Licenses'}
        placeHolderText={'ISA Document License'}
        options={ISA_DOCUMENT_LICENSES}
        multiple
      />
      <SelectBar
        label={'Contacts'}
        placeHolderText={'Contact'}
        options={CONTACTS}
        multiple
      />
    </Collapsible>
  );
};

export default InvestigationFormGeneral;
